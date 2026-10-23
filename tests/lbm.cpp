#include "../cpp/core.hpp"

#include <algorithm>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <iostream>
#include <stdexcept>

using lmfluid::LBM;
static constexpr double PI = 3.141592653589793238462643383279502884;
static int checks = 0;
static void check(bool condition, const char *what) {
  ++checks;
  if (!condition)
    throw std::runtime_error(what);
}
static double rel(double a, double b) {
  return std::abs(a - b) / std::max(1e-15, std::abs(b));
}
static bool throws_bad(void (*fn)()) {
  try {
    fn();
  } catch (const std::exception &) {
    return true;
  }
  return false;
}

static void basic_and_determinism() {
  auto a = LBM::periodic(16, 12, 0.08), b = LBM::periodic(16, 12, 0.08);
  a.init("0.03", "0.01").run(200);
  b.init("0.03", "0.01").run(200);
  check(rel(a.mass(), 192.0) < 1e-12, "periodic mass conservation");
  check(rel(a.velocity(4, 7, 0), 0.03) < 1e-12 &&
            rel(a.velocity(4, 7, 1), 0.01) < 1e-12,
        "periodic momentum");
  check(a.steps() == 200 && a.sample().data == b.sample().data,
        "independent runs deterministic");
  auto split = LBM::periodic(16, 12, 0.08);
  split.init("0.03", "0.01").run(100).run(100);
  check(split.sample().data == a.sample().data, "split and full run identical");
  auto s = a.sample();
  check(s.rows == 192 && s.cols == 6 && std::abs(s.get(0, 0) - .5) < 1e-14 &&
            std::abs(s.get(1, 1) - .5) < 1e-14,
        "sample schema and cell centers");
  a.init("0.01*t", "0");
  check(a.steps() == 0 && std::abs(a.velocity(3, 4, 0)) < 1e-14,
        "init resets time and state");
}

static void shear_decay() {
  constexpr int H = 64;
  constexpr double nu = .06, amp = .001;
  auto s = LBM::periodic(24, H, nu);
  s.init("0.001*sin(2*pi*y/64)", "0").run(100);
  double measured = 0,
         expected = amp * std::exp(-nu * std::pow(2 * PI / H, 2) * 100);
  for (int y = 0; y < H; ++y)
    measured += s.velocity(7, y, 0) * std::sin(2 * PI * (y + .5) / H) * 2.0 / H;
  check(rel(measured, expected) < .01, "periodic shear-wave viscous decay");
}

static void poiseuille() {
  constexpr int nx = 24, ny = 18;
  constexpr double nu = .1, g = 1e-5;
  auto c = LBM::channel(nx, ny, nu, g);
  c.init("0", "0").run(10000);
  double e2 = 0, r2 = 0;
  for (int y = 0; y < ny; ++y) {
    double exact = g * (y + .5) * (ny - (y + .5)) / (2 * nu);
    double got = c.velocity(nx / 2, y, 0);
    e2 += (got - exact) * (got - exact);
    r2 += exact * exact;
    check(std::abs(c.velocity(nx / 2, y, 1)) < 1e-10, "channel cross velocity");
  }
  check(std::sqrt(e2 / r2) < .01, "Poiseuille profile");
  check(rel(c.mass(), nx * ny) < 1e-9, "channel mass conservation");
}

// Ghia, Ghia & Shin, J. Comput. Phys. 48 (1982), Re=100, tables I and II: u on the vertical
// centerline (y, u/U) and v on the horizontal centerline (x, v/U).
static const double ghia_u[][2] = {
    {.9766, .84123},  {.9688, .78871},  {.9609, .73722},  {.9531, .68717},
    {.8516, .23151},  {.7344, .00332},  {.6172, -.13641}, {.5, -.20581},
    {.4531, -.21090}, {.2813, -.15662}, {.1719, -.10150}, {.1016, -.06434},
    {.0703, -.04775}, {.0625, -.04192}, {.0547, -.03717}};
static const double ghia_v[][2] = {
    {.9688, -.05906}, {.9609, -.07391}, {.9531, -.08864}, {.9453, -.10313},
    {.9063, -.16914}, {.8594, -.22445}, {.8047, -.24533}, {.5, .05454},
    {.2344, .17527},  {.2266, .17507},  {.1563, .16077},  {.0938, .12317},
    {.0781, .10890},  {.0703, .10091},  {.0625, .09233}};

// Linear interpolation along a centerline (average of the two middle columns/rows) between cell
// centers at (j + .5) / n; s is the normalized position along the line.
static double centerline(const LBM &c, int n, double s, int component) {
  double p = s * n - .5;
  int j = (int)std::floor(p);
  double a = p - j;
  auto at = [&](int k) {
    k = std::max(0, std::min(n - 1, k));
    return component == 0
               ? .5 * (c.velocity(n / 2 - 1, k, 0) + c.velocity(n / 2, k, 0))
               : .5 * (c.velocity(k, n / 2 - 1, 1) + c.velocity(k, n / 2, 1));
  };
  return (1 - a) * at(j) + a * at(j + 1);
}

static void ghia_cavity() {
  constexpr int n = 32;
  constexpr double U = .08;
  auto c = LBM::cavity(n, n, U * n / 100.0, U);
  auto t0 = std::chrono::steady_clock::now();
  c.init("0", "0").run(29000);
  double before = centerline(c, n, .5, 0);
  c.run(1000);
  double center = .25 * (c.velocity(15, 15, 0) + c.velocity(16, 15, 0) +
                         c.velocity(15, 16, 0) + c.velocity(16, 16, 0));
  double ratio = center / U;
  double sec =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - t0)
          .count();
  check(std::abs(centerline(c, n, .5, 0) - before) < 1e-6 * U,
        "cavity reached a steady state");
  check(rel(c.mass(), n * n) < 1e-9, "cavity mass conservation");
  check(std::abs(ratio + .20581) < .015, "Ghia Re=100 cavity center velocity");
  // Whole profiles: 30 tabulated points. Measured max deviation 0.0093 (u) / 0.0074 (v) at
  // 32x32, 0.0056 / 0.0034 at 64x64 (second-order approach to the benchmark).
  double du = 0, dv = 0;
  for (const auto &g : ghia_u)
    du = std::max(du, std::abs(centerline(c, n, g[0], 0) / U - g[1]));
  for (const auto &g : ghia_v)
    dv = std::max(dv, std::abs(centerline(c, n, g[0], 1) / U - g[1]));
  check(du < .015, "Ghia Re=100 u profile on the vertical centerline");
  check(dv < .015, "Ghia Re=100 v profile on the horizontal centerline");
  std::cout << "ghia center/U=" << ratio << ", max profile deviation u " << du
            << " v " << dv << " (30k steps, " << sec << " s)\n";
}

int main() {
  try {
    basic_and_determinism();
    shear_decay();
    poiseuille();
    ghia_cavity();
    check(throws_bad([] { (void)LBM::periodic(2, 4, .1); }),
          "small grid rejected");
    check(throws_bad([] { (void)LBM::periodic(4, 4, 0.0); }),
          "nonpositive viscosity rejected");
    check(throws_bad([] { (void)LBM::periodic(4, 4, DBL_MAX); }),
          "unrepresentable tau rejected");
    check(throws_bad([] { (void)LBM::cavity(4, 4, .1, .2); }),
          "fast lid rejected");
    check(throws_bad([] { (void)LBM::channel(4, 4, .1, .03); }),
          "large forcing rejected");
    check(throws_bad([] { (void)LBM::periodic(4, 4, 1e-300); }),
          "unresolved viscosity rejected");
    std::cout << "lbm checks: " << checks << " passed\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "lbm failure: " << e.what() << "\n";
    return 1;
  }
}
