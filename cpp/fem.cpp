#include "core.hpp"
#include "expr.hpp"

#include "mfem.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>

namespace lmfluid {

struct MeshData {
  std::unique_ptr<mfem::Mesh> mesh;
};

struct FieldData {
  std::shared_ptr<MeshData> owner;
  std::unique_ptr<mfem::FiniteElementCollection> fec;
  std::unique_ptr<mfem::FiniteElementSpace> fes;
  std::unique_ptr<mfem::GridFunction> gf;
  int vdim = 1;
  double current_time = 0.0;
};

namespace {

using mfem::real_t;

static void check(bool ok, const std::string &s) {
  if (!ok)
    throw std::runtime_error("lmfluid: " + s);
}

static std::vector<std::string> xyz_names() { return {"x", "y", "z", "t"}; }

struct ExprCoefficient : mfem::Coefficient {
  Expr expr;
  const double *time;
  bool positive;
  ExprCoefficient(const std::string &s, const double *tp = nullptr,
                  bool pos = false)
      : expr(s, xyz_names()), time(tp), positive(pos) {}
  real_t Eval(mfem::ElementTransformation &T,
              const mfem::IntegrationPoint &ip) override {
    mfem::Vector p(3);
    T.Transform(ip, p);
    double v[4] = {p(0), p.Size() > 1 ? p(1) : 0.0, p.Size() > 2 ? p(2) : 0.0,
                   time ? *time : 0.0};
    double result = expr.eval(v);
    if (!std::isfinite(result) || (positive && result <= 0.0))
      throw std::runtime_error(
          positive
              ? "lmfluid: diffusion coefficient must be finite and positive"
              : "lmfluid: coefficient evaluated to a nonfinite value");
    return result;
  }
};

struct PiecewiseCoefficient : mfem::Coefficient {
  std::map<int, std::unique_ptr<ExprCoefficient>> pieces;
  std::unique_ptr<ExprCoefficient> fallback;
  const double *time;
  PiecewiseCoefficient(const std::vector<Boundary> &bs, const double *tp)
      : time(tp) {
    for (const auto &b : bs) {
      if (b.attribute == 0)
        fallback = std::make_unique<ExprCoefficient>(b.x, tp);
      else
        pieces[b.attribute] = std::make_unique<ExprCoefficient>(b.x, tp);
    }
    if (!fallback)
      fallback = std::make_unique<ExprCoefficient>("0", tp);
  }
  real_t Eval(mfem::ElementTransformation &T,
              const mfem::IntegrationPoint &ip) override {
    auto it = pieces.find(T.Attribute);
    return (it == pieces.end() ? fallback : it->second)->Eval(T, ip);
  }
};

struct ExprVectorCoefficient : mfem::VectorCoefficient {
  std::vector<std::unique_ptr<ExprCoefficient>> exprs;
  ExprVectorCoefficient(int d, const std::vector<std::string> &s,
                        const double *tp = nullptr)
      : mfem::VectorCoefficient(d) {
    for (int i = 0; i < d; ++i)
      exprs.emplace_back(std::make_unique<ExprCoefficient>(
          i < (int)s.size() ? s[i] : "0", tp));
  }
  void Eval(mfem::Vector &v, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override {
    v.SetSize(vdim);
    for (int i = 0; i < vdim; ++i)
      v(i) = exprs[i]->Eval(T, ip);
  }
};

struct PiecewiseVectorCoefficient : mfem::VectorCoefficient {
  std::map<int, std::unique_ptr<ExprVectorCoefficient>> pieces;
  std::unique_ptr<ExprVectorCoefficient> fallback;
  PiecewiseVectorCoefficient(int dim, const std::vector<Boundary> &bs,
                             const double *tp)
      : mfem::VectorCoefficient(dim) {
    for (const auto &b : bs) {
      std::vector<std::string> s = {b.x, b.y, b.z};
      if (b.attribute == 0)
        fallback = std::make_unique<ExprVectorCoefficient>(dim, s, tp);
      else
        pieces[b.attribute] =
            std::make_unique<ExprVectorCoefficient>(dim, s, tp);
    }
    if (!fallback)
      fallback = std::make_unique<ExprVectorCoefficient>(
          dim, std::vector<std::string>{"0", "0", "0"}, tp);
  }
  void Eval(mfem::Vector &v, mfem::ElementTransformation &T,
            const mfem::IntegrationPoint &ip) override {
    auto it = pieces.find(T.Attribute);
    (it == pieces.end() ? fallback : it->second)->Eval(v, T, ip);
  }
};

static std::shared_ptr<MeshData> mesh_owner(const Mesh &m) {
  check(m.data && m.data->mesh, "mesh is empty");
  return m.data;
}

static mfem::Array<int> boundary_marker(const mfem::Mesh &m,
                                        const std::vector<Boundary> &bs) {
  check(m.bdr_attributes.Size() > 0, "mesh has no boundary attributes");
  mfem::Array<int> mark(m.bdr_attributes.Max());
  mark = 0;
  for (const auto &b : bs) {
    if (b.attribute == 0)
      mark = 1;
    else {
      check(b.attribute > 0 && b.attribute <= mark.Size() &&
                m.bdr_attributes.Find(b.attribute) >= 0,
            "boundary attribute does not exist");
      mark[b.attribute - 1] = 1;
    }
  }
  return mark;
}

static std::unique_ptr<FieldData>
make_field(const std::shared_ptr<MeshData> &mo, int order, int vdim) {
  int dim = mo->mesh->Dimension();
  check(order >= 1 && order <= 3, "invalid FEM order");
  auto f = std::make_unique<FieldData>();
  f->owner = mo;
  f->vdim = vdim;
  f->fec = std::make_unique<mfem::H1_FECollection>(order, dim);
  f->fes = std::make_unique<mfem::FiniteElementSpace>(
      mo->mesh.get(), f->fec.get(), vdim, mfem::Ordering::byVDIM);
  check(f->fes->GetVSize() > 0 && f->fes->GetVSize() <= 1000000,
        "finite element space is too large");
  f->gf = std::make_unique<mfem::GridFunction>(f->fes.get());
  *f->gf = 0.0;
  return f;
}

static void solve_linear(mfem::BilinearForm &a, mfem::LinearForm &b,
                         mfem::GridFunction &x, const mfem::Array<int> &ess,
                         double rtol, int maxit) {
  mfem::OperatorPtr A;
  mfem::Vector B, X;
  check(std::isfinite(rtol) && rtol > 0 && rtol < 1 && maxit > 0 &&
            maxit <= 100000,
        "invalid solver tolerance or iteration limit");
  a.FormLinearSystem(ess, x, b, A, X, B);
  auto *sp = A.As<mfem::SparseMatrix>();
  check(sp != nullptr,
        "serial FEM system was not assembled as a sparse matrix");
  mfem::GSSmoother M(*sp);
  mfem::CGSolver cg;
  cg.SetPreconditioner(M);
  cg.SetOperator(*A);
  cg.SetMaxIter(maxit);
  cg.SetRelTol(rtol);
  cg.SetAbsTol(0.0);
  cg.SetPrintLevel(0);
  cg.Mult(B, X);
  check(cg.GetConverged(), "linear solve did not converge (residual " +
                               std::to_string(cg.GetFinalNorm()) + ")");
  a.RecoverFEMSolution(X, B, x);
  for (int i = 0; i < x.Size(); ++i)
    check(std::isfinite(x[i]), "linear solve returned a nonfinite field");
}

} // namespace

Mesh Mesh::rect(int nx, int ny, double lx, double ly) {
  check(nx > 0 && ny > 0 && std::isfinite(lx) && std::isfinite(ly) && lx > 0 &&
            ly > 0,
        "invalid rectangular mesh");
  check(static_cast<size_t>(nx) * ny <= 100000, "mesh element limit exceeded");
  Mesh out;
  out.data = std::make_shared<MeshData>();
  out.data->mesh = std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian2D(
      nx, ny, mfem::Element::QUADRILATERAL, false, lx, ly, false));
  out.data->mesh->SetAttributes();
  return out;
}

Mesh Mesh::box(int nx, int ny, int nz, double lx, double ly, double lz) {
  check(nx > 0 && ny > 0 && nz > 0 && lx > 0 && ly > 0 && lz > 0 &&
            std::isfinite(lx) && std::isfinite(ly) && std::isfinite(lz),
        "invalid box mesh");
  check(static_cast<size_t>(nx) * ny <= 100000 / static_cast<size_t>(nz),
        "mesh element limit exceeded");
  Mesh out;
  out.data = std::make_shared<MeshData>();
  out.data->mesh = std::make_unique<mfem::Mesh>(mfem::Mesh::MakeCartesian3D(
      nx, ny, nz, mfem::Element::HEXAHEDRON, lx, ly, lz, false));
  out.data->mesh->SetAttributes();
  return out;
}

Mesh Mesh::load(const std::string &path) {
  check(!path.empty(), "mesh path is empty");
  std::ifstream input(path);
  check((bool)input, "cannot open mesh: " + path);
  Mesh out;
  out.data = std::make_shared<MeshData>();
  out.data->mesh.reset(new mfem::Mesh(path.c_str(), 1, 1, true));
  check(out.data->mesh->Dimension() == 2 || out.data->mesh->Dimension() == 3,
        "only 2D and 3D meshes are supported");
  check(out.data->mesh->Dimension() == out.data->mesh->SpaceDimension(),
        "embedded surface meshes are unsupported");
  check(out.data->mesh->GetNE() > 0 && out.data->mesh->GetNE() <= 100000,
        "mesh element limit exceeded");
  return out;
}

int Mesh::dimension() const {
  return data && data->mesh ? data->mesh->Dimension() : 0;
}
int Mesh::elements() const {
  return data && data->mesh ? data->mesh->GetNE() : 0;
}
void Mesh::save(const std::string &path) const {
  check(data && data->mesh, "mesh is empty");
  std::ofstream os(path);
  check((bool)os, "cannot open mesh output: " + path);
  os << std::setprecision(17);
  data->mesh->Print(os);
  os.flush();
  check((bool)os, "cannot write mesh output: " + path);
}

double Field::value(double x, double y, double z, int component) const {
  check(data && data->gf, "field is empty");
  check(component >= 0 && component < data->vdim,
        "field component out of range");
  int dim = data->owner->mesh->Dimension();
  check(std::isfinite(x) && std::isfinite(y) && std::isfinite(z),
        "point is nonfinite");
  if (dim == 2)
    check(std::abs(z) <= 1e-12, "2D point must have z=0");
  mfem::DenseMatrix p(dim, 1);
  p(0, 0) = x;
  if (dim > 1)
    p(1, 0) = y;
  if (dim > 2)
    p(2, 0) = z;
  mfem::Array<int> elems;
  mfem::Array<mfem::IntegrationPoint> ips;
  data->owner->mesh->FindPoints(p, elems, ips, false);
  check(elems.Size() == 1 && elems[0] >= 0, "point is outside mesh");
  auto *T = data->owner->mesh->GetElementTransformation(elems[0]);
  if (data->vdim == 1)
    return data->gf->GetValue(*T, ips[0]);
  mfem::Vector v;
  data->gf->GetVectorValue(*T, ips[0], v);
  return v(component);
}

double Field::l2_error(const std::string &s) const {
  check(data && data->gf && data->vdim == 1,
        "scalar L2 error requires a scalar field");
  ExprCoefficient c(s, &data->current_time);
  return data->gf->ComputeL2Error(c);
}

double Field::l2_error_vector(const std::string &x, const std::string &y,
                              const std::string &z) const {
  check(data && data->gf && data->vdim > 1,
        "vector L2 error requires a vector field");
  int d = data->owner->mesh->Dimension();
  ExprVectorCoefficient c(d, {x, y, z}, &data->current_time);
  return data->gf->ComputeL2Error(c);
}

Matrix Field::sample(int nx, int ny, int nz) const {
  check(data && data->gf, "field is empty");
  int d = data->owner->mesh->Dimension();
  check(nx >= 2 && nx <= 4096 && ny >= 2 && ny <= 4096 && nz <= 4096 &&
            (d == 2 ? nz == 1 : nz >= 2),
        "invalid sample dimensions");
  size_t count = static_cast<size_t>(nx) * ny * (d == 3 ? nz : 1);
  check(count <= 5000000 &&
            count <= std::numeric_limits<size_t>::max() / (3 + data->vdim),
        "sample is too large");
  mfem::Vector lo, hi;
  data->owner->mesh->GetBoundingBox(lo, hi);
  int zcount = d == 3 ? nz : 1;
  Matrix out;
  out.rows = nx * ny * zcount;
  out.cols = 3 + data->vdim;
  out.data.resize((size_t)out.rows * out.cols);
  int row = 0;
  for (int k = 0; k < zcount; ++k)
    for (int j = 0; j < ny; ++j)
      for (int i = 0; i < nx; ++i, ++row) {
        double x = nx == 1 ? (lo(0) + hi(0)) / 2
                           : lo(0) + (hi(0) - lo(0)) * i / (nx - 1.0);
        double y = ny == 1 ? (lo(1) + hi(1)) / 2
                           : lo(1) + (hi(1) - lo(1)) * j / (ny - 1.0);
        double z = d == 3 ? (nz == 1 ? (lo(2) + hi(2)) / 2
                                     : lo(2) + (hi(2) - lo(2)) * k / (nz - 1.0))
                          : 0;
        out.data[row * out.cols] = x;
        out.data[row * out.cols + 1] = y;
        out.data[row * out.cols + 2] = z;
        for (int c = 0; c < data->vdim; ++c)
          out.data[row * out.cols + 3 + c] = value(x, y, z, c);
      }
  return out;
}

int Field::dofs() const {
  return data && data->fes ? data->fes->GetVSize() : 0;
}
int Field::components() const { return data ? data->vdim : 0; }
double Field::time() const { return data ? data->current_time : 0.0; }

Problem &Problem::diffusion(const std::string &s) {
  coefficient = s;
  return *this;
}
Problem &Problem::source(const std::string &s) {
  source_expr = s;
  return *this;
}
Problem &Problem::initial(const std::string &s) {
  initial_expr = s;
  return *this;
}
Problem &Problem::dirichlet(int a, const std::string &s) {
  if (a == 0)
    boundaries.clear();
  else
    boundaries.erase(
        std::remove_if(boundaries.begin(), boundaries.end(),
                       [](const Boundary &b) { return b.attribute == 0; }),
        boundaries.end());
  boundaries.push_back({a, s, "0", "0"});
  return *this;
}
Problem &Problem::displacement(int a, const std::string &x,
                               const std::string &y, const std::string &z) {
  if (a == 0)
    boundaries.clear();
  else
    boundaries.erase(
        std::remove_if(boundaries.begin(), boundaries.end(),
                       [](const Boundary &b) { return b.attribute == 0; }),
        boundaries.end());
  boundaries.push_back({a, x, y, z});
  return *this;
}
Problem &Problem::force(const std::string &x, const std::string &y,
                        const std::string &z) {
  force_x = x;
  force_y = y;
  force_z = z;
  return *this;
}
Problem &Problem::lame(double l, double m) {
  check(m > 0 && std::isfinite(l) && std::isfinite(m),
        "invalid Lame parameters");
  lambda = l;
  mu = m;
  return *this;
}
Problem &Problem::tol(double r, int n) {
  check(r > 0 && r < 1 && n > 0, "invalid solver tolerance or iteration limit");
  relative_tolerance = r;
  max_iterations = n;
  return *this;
}

Field Problem::solve() const {
  check(mesh.data && mesh.data->mesh, "problem mesh is empty");
  int dim = mesh.dimension();
  check(order >= 1 && order <= 3, "finite element order must be 1..3");
  check(kind == "poisson" || kind == "elasticity",
        "solve supports poisson or elasticity");
  if (kind == "elasticity")
    check(mu > 0 && lambda + 2.0 * mu / dim > 0 && std::isfinite(lambda) &&
              std::isfinite(mu),
          "invalid Lame parameters");
  auto mo = mesh_owner(mesh);
  int vd = kind == "elasticity" ? dim : 1;
  auto fd = make_field(mo, order, vd);
  mfem::BilinearForm a(fd->fes.get());
  mfem::LinearForm b(fd->fes.get());
  ExprCoefficient k(coefficient, nullptr, true), src(source_expr);
  std::unique_ptr<mfem::ConstantCoefficient> lc, mc;
  std::unique_ptr<ExprVectorCoefficient> force_coeff;
  if (kind == "poisson") {
    a.AddDomainIntegrator(new mfem::DiffusionIntegrator(k));
    b.AddDomainIntegrator(new mfem::DomainLFIntegrator(src));
  } else {
    lc = std::make_unique<mfem::ConstantCoefficient>(lambda);
    mc = std::make_unique<mfem::ConstantCoefficient>(mu);
    force_coeff = std::make_unique<ExprVectorCoefficient>(
        dim, std::vector<std::string>{force_x, force_y, force_z});
    a.AddDomainIntegrator(new mfem::ElasticityIntegrator(*lc, *mc));
    b.AddDomainIntegrator(new mfem::VectorDomainLFIntegrator(*force_coeff));
  }
  a.Assemble();
  b.Assemble();
  auto marker = boundary_marker(*mo->mesh, boundaries);
  mfem::Array<int> ess;
  fd->fes->GetEssentialTrueDofs(marker, ess);
  check(ess.Size() > 0, "steady FEM requires Dirichlet boundary data");
  if (kind == "poisson") {
    PiecewiseCoefficient bc(boundaries, &fd->current_time);
    fd->gf->ProjectBdrCoefficient(bc, marker);
  } else {
    PiecewiseVectorCoefficient bc(dim, boundaries, &fd->current_time);
    fd->gf->ProjectBdrCoefficient(bc, marker);
  }
  solve_linear(a, b, *fd->gf, ess, relative_tolerance, max_iterations);
  Field out;
  out.data = std::move(fd);
  return out;
}

Field Problem::evolve(double tend, double dt) const {
  check(kind == "heat", "evolve supports heat only");
  check(tend >= 0 && dt > 0 && std::isfinite(tend) && std::isfinite(dt),
        "invalid heat time interval or step");
  check(mesh.data && mesh.data->mesh, "problem mesh is empty");
  check(order >= 1 && order <= 3, "finite element order must be 1..3");
  auto mo = mesh_owner(mesh);
  auto fd = make_field(mo, order, 1);
  ExprCoefficient k(coefficient, &fd->current_time, true),
      src(source_expr, &fd->current_time),
      init(initial_expr, &fd->current_time);
  fd->gf->ProjectCoefficient(init);
  auto mark = boundary_marker(*mo->mesh, boundaries);
  mfem::Array<int> ess;
  fd->fes->GetEssentialTrueDofs(mark, ess);
  double t = 0;
  double nreal = std::ceil(tend / dt);
  check(std::isfinite(nreal) && nreal <= 1000000, "too many time steps");
  int n = (int)nreal;
  PiecewiseCoefficient initial_bc(boundaries, &fd->current_time);
  fd->gf->ProjectBdrCoefficient(initial_bc, mark);
  for (int step = 0; step < n && t < tend; ++step) {
    const double target =
        step + 1 == n ? tend : std::min(tend, (step + 1) * dt);
    double h = target - t;
    check(h > 0 && std::isfinite(h) && std::isfinite(1.0 / h) && t + h > t,
          "time stepping made no progress");
    t = target;
    fd->current_time = t;
    mfem::BilinearForm a(fd->fes.get());
    mfem::ConstantCoefficient invh(1.0 / h);
    a.AddDomainIntegrator(new mfem::DiffusionIntegrator(k));
    a.AddDomainIntegrator(new mfem::MassIntegrator(invh));
    mfem::LinearForm b(fd->fes.get());
    b.AddDomainIntegrator(new mfem::DomainLFIntegrator(src));
    mfem::GridFunction old(*fd->gf);
    mfem::GridFunctionCoefficient oldc(&old);
    mfem::ProductCoefficient old_over_h(oldc, invh);
    b.AddDomainIntegrator(new mfem::DomainLFIntegrator(old_over_h));
    a.Assemble();
    b.Assemble();
    PiecewiseCoefficient bc(boundaries, &fd->current_time);
    fd->gf->ProjectBdrCoefficient(bc, mark);
    solve_linear(a, b, *fd->gf, ess, relative_tolerance, max_iterations);
  }
  Field out;
  out.data = std::move(fd);
  return out;
}

} // namespace lmfluid
