#include "core.hpp"
#include "expr.hpp"

#include <algorithm>
#include <cfloat>
#include <climits>
#include <cmath>
#include <iostream>
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
  // solid[c] != 0 marks an obstacle cell: no collision, halfway bounceback on
  // every fluid link that points into it, zero velocity and unit density.
  std::vector<char> solid;
  bool any_solid = false, warned_mach = false;
  // momentum exchanged with obstacle cells during the most recent step, and its
  // mean over the most recent run() call (which filters acoustic noise)
  double solid_fx = 0, solid_fy = 0, mean_fx = 0, mean_fy = 0;
  int nsteps = 0;

  std::size_t cell(int x, int y) const { return (std::size_t)y * nx + x; }
  std::size_t at(int x, int y, int q) const { return cell(x, y) * Q + q; }
  void validate_state() const {
    for (std::size_t k = 0; k < rho.size(); ++k) {
      if (solid[k])
        continue;
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
        if (solid[c]) {
          rho[c] = 1;
          ux[c] = uy[c] = 0;
          continue;
        }
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
  // rho/ux/uy are current on entry: construction, init(), obstacle() and the
  // end of every step() refresh them.
  void step() {
    const double tau = 0.5 + 3.0 * nu, omega = 1.0 / tau;
    std::fill(next.begin(), next.end(), 0.0);
    double fx = 0, fy = 0;
    for (int y = 0; y < ny; ++y)
      for (int x = 0; x < nx; ++x) {
        std::size_t c = cell(x, y);
        if (solid[c])
          continue;
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
          if (outside(xx, yy))
            reflected(x, y, q, post);
          else if (solid[cell(xx, yy)]) {
            // Stationary halfway bounceback; the obstacle receives 2*f*e.
            next[at(x, y, opp[q])] += post;
            fx += 2.0 * post * ex[q];
            fy += 2.0 * post * ey[q];
          } else
            next[at(xx, yy, q)] += post;
        }
      }
    // Keep solid populations at rest so their (unused) values stay finite.
    if (any_solid)
      for (std::size_t c = 0; c < solid.size(); ++c)
        if (solid[c])
          for (int q = 0; q < Q; ++q)
            next[c * Q + q] = wt[q];
    solid_fx = fx;
    solid_fy = fy;
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
  if (0.5 + 3.0 * nu < 0.51)
    std::cerr << "lmfluid: warning: LBM relaxation time tau=" << 0.5 + 3.0 * nu
              << " is close to 0.5; BGK is fragile there (raise the viscosity "
                 "or refine the grid)\n";
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
  d->solid.assign(n, 0);
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
  data->mean_fx = data->mean_fy = 0;
  Expr exx(sx, {"x", "y", "t"}), exy(sy, {"x", "y", "t"});
  double v[3];
  for (int y = 0; y < data->ny; ++y)
    for (int x = 0; x < data->nx; ++x) {
      std::size_t c = data->cell(x, y);
      if (data->solid[c])
        continue;
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
  double sx = 0, sy = 0;
  for (int i = 0; i < steps; ++i) {
    data->step();
    sx += data->solid_fx;
    sy += data->solid_fy;
  }
  if (steps > 0) {
    data->mean_fx = sx / steps;
    data->mean_fy = sy / steps;
  }
  if (!data->warned_mach && mach() > 0.3) {
    data->warned_mach = true;
    std::cerr << "lmfluid: warning: LBM Mach number reached " << mach()
              << " (above 0.3, speed above 0.17 lattice units); "
                 "compressibility errors grow, so lower the speeds\n";
  }
  return *this;
}
LBM &LBM::obstacle(const std::string &mask) {
  if (!data)
    throw std::logic_error("LBM object is empty");
  Expr e(mask, {"x", "y"});
  double v[2];
  std::size_t fluid = 0;
  std::vector<char> next = data->solid;
  for (int y = 0; y < data->ny; ++y)
    for (int x = 0; x < data->nx; ++x) {
      v[0] = x + 0.5;
      v[1] = y + 0.5;
      double m = e.eval(v);
      if (!std::isfinite(m))
        throw std::invalid_argument("LBM obstacle expression is not finite");
      std::size_t c = data->cell(x, y);
      if (m > 0)
        next[c] = 1;
      fluid += !next[c];
    }
  if (fluid == 0)
    throw std::invalid_argument("LBM obstacle leaves no fluid cells");
  data->solid.swap(next);
  data->any_solid = std::find(data->solid.begin(), data->solid.end(), 1) !=
                    data->solid.end();
  for (std::size_t c = 0; c < data->solid.size(); ++c)
    if (data->solid[c])
      for (int q = 0; q < Q; ++q)
        data->f[c * Q + q] = wt[q];
  data->macros();
  return *this;
}
bool LBM::is_solid(int ix, int iy) const {
  if (!data)
    throw std::logic_error("LBM object is empty");
  if (ix < 0 || ix >= data->nx || iy < 0 || iy >= data->ny)
    throw std::out_of_range("LBM solid index out of range");
  return data->solid[data->cell(ix, iy)] != 0;
}
int LBM::solid_cells() const {
  if (!data)
    throw std::logic_error("LBM object is empty");
  return (int)std::count(data->solid.begin(), data->solid.end(), 1);
}
double LBM::solid_force(int component) const {
  if (!data)
    throw std::logic_error("LBM object is empty");
  if (component < 0 || component > 1)
    throw std::out_of_range("LBM force component must be 0 or 1");
  return component ? data->mean_fy : data->mean_fx;
}
double LBM::tau() const {
  if (!data)
    throw std::logic_error("LBM object is empty");
  return 0.5 + 3.0 * data->nu;
}
double LBM::mach() const { return max_speed() * std::sqrt(3.0); }
Matrix LBM::vorticity() const {
  if (!data)
    throw std::logic_error("LBM object is empty");
  const LbmData &d = *data;
  // Central differences; a wall or obstacle neighbour contributes its no-slip
  // zero velocity at half a cell, so the stencil stays one-sided there.
  auto vel = [&](int x, int y, int comp, bool &wall) -> double {
    wall = false;
    if (d.x_periodic)
      x = (x + d.nx) % d.nx;
    if (d.y_periodic)
      y = (y + d.ny) % d.ny;
    if (d.outside(x, y)) {
      wall = true;
      // only the moving lid has a nonzero wall velocity
      return (d.moving_top && y >= d.ny && comp == 0) ? d.lid : 0.0;
    }
    std::size_t c = d.cell(x, y);
    if (d.solid[c]) {
      wall = true;
      return 0.0;
    }
    return comp ? d.uy[c] : d.ux[c];
  };
  auto deriv = [&](int x, int y, int dx, int dy, int comp) {
    bool wp, wm, w0;
    double up = vel(x + dx, y + dy, comp, wp), um = vel(x - dx, y - dy, comp, wm);
    double u0 = vel(x, y, comp, w0);
    double hp = wp ? 0.5 : 1.0, hm = wm ? 0.5 : 1.0;
    return (up - u0) / hp * hm / (hp + hm) + (u0 - um) / hm * hp / (hp + hm);
  };
  Matrix m;
  m.rows = d.nx * d.ny;
  m.cols = 4;
  m.data.reserve((std::size_t)m.rows * m.cols);
  for (int y = 0; y < d.ny; ++y)
    for (int x = 0; x < d.nx; ++x) {
      std::size_t c = d.cell(x, y);
      double w = d.solid[c] ? 0.0 : deriv(x, y, 1, 0, 1) - deriv(x, y, 0, 1, 0);
      m.data.insert(m.data.end(),
                    {(double)x + 0.5, (double)y + 0.5, w, d.solid[c] ? 1.0 : 0.0});
    }
  return m;
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
  double m = 0;
  for (std::size_t i = 0; i < data->rho.size(); ++i)
    if (!data->solid[i])
      m += data->rho[i];
  return m;
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
