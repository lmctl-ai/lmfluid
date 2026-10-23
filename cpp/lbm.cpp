#include "core.hpp"
#include "expr.hpp"

#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>

namespace lmfluid {
namespace {
constexpr int Q = 9;
constexpr int ex[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
constexpr int ey[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
constexpr int opp[Q] = {0, 3, 4, 1, 2, 7, 8, 5, 6};
constexpr double wt[Q] = {4.0 / 9.0,  1.0 / 9.0,  1.0 / 9.0,
                          1.0 / 9.0,  1.0 / 9.0,  1.0 / 36.0,
                          1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0};

} // namespace

struct LbmData {
  int nx = 0, ny = 0;
  double nu = 0, force_x = 0, lid = 0;
  bool x_periodic = false, y_periodic = false, moving_top = false;
  std::vector<double> f, next, rho, ux, uy;
  int nsteps = 0;

  std::size_t cell(int x, int y) const { return (std::size_t)y * nx + x; }
  std::size_t at(int x, int y, int q) const { return cell(x, y) * Q + q; }
  void validate_state() const {
    for (std::size_t k = 0; k < rho.size(); ++k) {
      if (!(std::isfinite(rho[k]) && rho[k] > 0 && std::isfinite(ux[k]) &&
            std::isfinite(uy[k])))
        throw std::runtime_error("LBM became numerically unstable (non-finite "
                                 "or non-positive state)");
      double u = std::hypot(ux[k], uy[k]);
      if (!(std::isfinite(u) && u < 0.45))
        throw std::runtime_error(
            "LBM became numerically unstable (Mach number too high)");
    }
  }
  void macros() {
    for (int y = 0; y < ny; ++y)
      for (int x = 0; x < nx; ++x) {
        std::size_t c = cell(x, y);
        double d = 0, px = 0, py = 0;
        for (int q = 0; q < Q; ++q) {
          double z = f[at(x, y, q)];
          d += z;
          px += z * ex[q];
          py += z * ey[q];
        }
        rho[c] = d;
        ux[c] = px / d + 0.5 * force_x;
        uy[c] = py / d;
      }
    validate_state();
  }
  bool outside(int x, int y) const {
    return x < 0 || x >= nx || y < 0 || y >= ny;
  }
  void reflected(int x, int y, int q, double val) {
    // Halfway bounceback. For a moving top wall, the standard momentum
    // correction is -6*w*rho*(e dot Uwall); see Zou & He (1997).
    double v = val;
    if (moving_top && ey[q] > 0 && y == ny - 1)
      v -= 6.0 * wt[q] * rho[cell(x, y)] * ex[q] * lid;
    next[at(x, y, opp[q])] += v;
  }
  void step() {
    macros();
    const double tau = 0.5 + 3.0 * nu, omega = 1.0 / tau;
    std::fill(next.begin(), next.end(), 0.0);
    for (int y = 0; y < ny; ++y)
      for (int x = 0; x < nx; ++x) {
        std::size_t c = cell(x, y);
        double d = rho[c], vx = ux[c], vy = uy[c];
        for (int q = 0; q < Q; ++q) {
          double eu = ex[q] * vx + ey[q] * vy;
          double feq =
              wt[q] * d *
              (1.0 + 3.0 * eu + 4.5 * eu * eu - 1.5 * (vx * vx + vy * vy));
          // Guo forcing uses F=rho*g; this keeps acceleration independent of
          // small compressibility fluctuations in rho.
          double F = d * force_x;
          double ef = ex[q] * F;
          double guo = wt[q] * (1.0 - 0.5 * omega) *
                       (3.0 * (ex[q] - vx) * F + 9.0 * eu * ef);
          double post = f[at(x, y, q)] - omega * (f[at(x, y, q)] - feq) + guo;
          int xx = x + ex[q], yy = y + ey[q];
          if (x_periodic) {
            if (xx < 0)
              xx += nx;
            if (xx >= nx)
              xx -= nx;
          }
          if (y_periodic) {
            if (yy < 0)
              yy += ny;
            if (yy >= ny)
              yy -= ny;
          }
          if (!outside(xx, yy))
            next[at(xx, yy, q)] += post;
          else
            reflected(x, y, q, post);
        }
      }
    f.swap(next);
    ++nsteps;
    macros();
  }
};

namespace {
void check_grid(int nx, int ny) {
  if (nx < 3 || ny < 3)
    throw std::invalid_argument("LBM grid must be at least 3 by 3");
  if ((std::size_t)nx * (std::size_t)ny > 1000000u)
    throw std::invalid_argument("LBM grid exceeds 1 million cells");
}
std::shared_ptr<LbmData> make(int nx, int ny, double nu) {
  check_grid(nx, ny);
  if (!std::isfinite(nu) || nu <= 0 || !std::isfinite(0.5 + 3.0 * nu) ||
      0.5 + 3.0 * nu <= 0.5)
    throw std::invalid_argument(
        "LBM viscosity must be finite, positive, and representable in tau");
  auto d = std::make_shared<LbmData>();
  d->nx = nx;
  d->ny = ny;
  d->nu = nu;
  std::size_t n = (std::size_t)nx * ny;
  d->f.assign(n * Q, 0.0);
  d->next.resize(n * Q);
  d->rho.assign(n, 1.0);
  d->ux.assign(n, 0.0);
  d->uy.assign(n, 0.0);
  for (std::size_t c = 0; c < n; ++c)
    for (int q = 0; q < Q; ++q)
      d->f[c * Q + q] = wt[q];
  return d;
}
void check_force(double g) {
  if (!std::isfinite(g) || std::abs(g) > 0.02)
    throw std::invalid_argument(
        "LBM acceleration must be finite and within +/-0.02");
}
void check_lid(double u) {
  if (!std::isfinite(u) || std::abs(u) >= 0.2)
    throw std::invalid_argument(
        "LBM lid speed must be finite and have magnitude below 0.2");
}
} // namespace

LBM LBM::channel(int nx, int ny, double viscosity, double acceleration) {
  check_force(acceleration);
  auto d = make(nx, ny, viscosity);
  d->x_periodic = true;
  d->force_x = acceleration;
  LBM r;
  r.data = std::move(d);
  r.init("0", "0");
  return r;
}
LBM LBM::periodic(int nx, int ny, double viscosity) {
  auto d = make(nx, ny, viscosity);
  d->x_periodic = d->y_periodic = true;
  LBM r;
  r.data = std::move(d);
  return r;
}
LBM LBM::cavity(int nx, int ny, double viscosity, double lid_speed) {
  check_lid(lid_speed);
  auto d = make(nx, ny, viscosity);
  d->moving_top = true;
  d->lid = lid_speed;
  LBM r;
  r.data = std::move(d);
  return r;
}
LBM &LBM::init(const std::string &sx, const std::string &sy) {
  if (!data)
    throw std::logic_error("LBM object is empty");
  data->nsteps = 0;
  Expr exx(sx, {"x", "y", "t"}), exy(sy, {"x", "y", "t"});
  double v[3];
  for (int y = 0; y < data->ny; ++y)
    for (int x = 0; x < data->nx; ++x) {
      std::size_t c = data->cell(x, y);
      v[0] = x + 0.5;
      v[1] = y + 0.5;
      v[2] = data->nsteps;
      double u = exx.eval(v), w = exy.eval(v);
      if (!std::isfinite(u) || !std::isfinite(w) || std::hypot(u, w) >= 0.35)
        throw std::invalid_argument(
            "LBM initial velocity is non-finite or exceeds stable Mach bound");
      data->rho[c] = 1;
      data->ux[c] = u;
      data->uy[c] = w;
      for (int q = 0; q < Q; ++q) {
        double eu = ex[q] * u + ey[q] * w;
        // The macroscopic velocity includes Guo's half-step force. Remove
        // that half impulse from the initialized populations so the requested
        // expression is the velocity reported at t=0.
        data->f[data->at(x, y, q)] =
            wt[q] * (1 + 3 * eu + 4.5 * eu * eu - 1.5 * (u * u + w * w)) -
            1.5 * wt[q] * ex[q] * data->force_x;
      }
    }
  data->macros();
  return *this;
}
LBM &LBM::run(int steps) {
  if (!data)
    throw std::logic_error("LBM object is empty");
  if (steps < 0 || steps > 10000000)
    throw std::invalid_argument("LBM steps must be between 0 and 10000000");
  if (steps > INT_MAX - data->nsteps)
    throw std::invalid_argument("LBM step count would overflow");
  if ((long double)steps * data->nx * data->ny > 500000000.0L)
    throw std::invalid_argument("LBM run exceeds work limit");
  for (int i = 0; i < steps; ++i)
    data->step();
  return *this;
}
double LBM::velocity(int ix, int iy, int component) const {
  if (!data)
    throw std::logic_error("LBM object is empty");
  if (ix < 0 || ix >= data->nx || iy < 0 || iy >= data->ny)
    throw std::out_of_range("LBM velocity index out of range");
  if (component < 0 || component > 1)
    throw std::out_of_range("LBM velocity component must be 0 or 1");
  return component ? data->uy[data->cell(ix, iy)]
                   : data->ux[data->cell(ix, iy)];
}
double LBM::density(int ix, int iy) const {
  if (!data)
    throw std::logic_error("LBM object is empty");
  if (ix < 0 || ix >= data->nx || iy < 0 || iy >= data->ny)
    throw std::out_of_range("LBM density index out of range");
  return data->rho[data->cell(ix, iy)];
}
double LBM::mass() const {
  if (!data)
    throw std::logic_error("LBM object is empty");
  return std::accumulate(data->rho.begin(), data->rho.end(), 0.0);
}
double LBM::max_speed() const {
  if (!data)
    throw std::logic_error("LBM object is empty");
  double m = 0;
  for (std::size_t i = 0; i < data->rho.size(); ++i)
    m = std::max(m, std::hypot(data->ux[i], data->uy[i]));
  return m;
}
double LBM::kinetic_energy() const {
  if (!data)
    throw std::logic_error("LBM object is empty");
  double e = 0;
  for (std::size_t i = 0; i < data->rho.size(); ++i)
    e += 0.5 * data->rho[i] *
         (data->ux[i] * data->ux[i] + data->uy[i] * data->uy[i]);
  return e;
}
int LBM::steps() const {
  if (!data)
    throw std::logic_error("LBM object is empty");
  return data->nsteps;
}
Matrix LBM::sample() const {
  if (!data)
    throw std::logic_error("LBM object is empty");
  Matrix m;
  m.rows = data->nx * data->ny;
  m.cols = 6;
  m.data.reserve((std::size_t)m.rows * m.cols);
  for (int y = 0; y < data->ny; ++y)
    for (int x = 0; x < data->nx; ++x) {
      auto c = data->cell(x, y);
      m.data.insert(m.data.end(), {(double)x + 0.5, (double)y + 0.5,
                                   data->rho[c], data->ux[c], data->uy[c],
                                   std::hypot(data->ux[c], data->uy[c])});
    }
  return m;
}
} // namespace lmfluid
