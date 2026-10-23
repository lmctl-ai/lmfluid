#include "../cpp/core.hpp"

#include <cmath>
#include <cstdio>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

using lmfluid::Field;
using lmfluid::Mesh;
using lmfluid::Problem;

static int checks = 0;
static void ok(bool condition, const char *what) {
  ++checks;
  if (!condition)
    throw std::runtime_error(what);
}
template <class F> static bool throws_bad(F &&fn) {
  try {
    fn();
  } catch (const std::exception &) {
    return true;
  }
  return false;
}

static double poisson_error(int n, int order) {
  Problem p;
  p.mesh = Mesh::rect(n, n, 1.0, 1.0);
  p.order = order;
  p.source("2*pi^2*sin(pi*x)*sin(pi*y)");
  return p.solve().l2_error("sin(pi*x)*sin(pi*y)");
}

static void convergence() {
  double p1a = poisson_error(4, 1), p1b = poisson_error(8, 1);
  double p2a = poisson_error(2, 2), p2b = poisson_error(4, 2);
  std::cout << "fem convergence errors: P1 " << p1a << " -> " << p1b
            << " (rate " << std::log(p1a / p1b) / std::log(2.0) << ")"
            << ", P2 " << p2a << " -> " << p2b << " (rate "
            << std::log(p2a / p2b) / std::log(2.0) << ")\n";
  ok(p1b < p1a / 3.0, "P1 L2 convergence order");
  ok(p2b < p2a / 6.0, "P2 L2 convergence order");
  ok(p1b < 0.02 && p2b < 0.002, "manufactured sine accuracy");
}

static void scalar_polynomial_and_coefficient() {
  Problem p;
  p.mesh = Mesh::rect(4, 4, 1.0, 1.0);
  p.order = 2;
  p.diffusion("2").source("-8").dirichlet(0, "x^2+y^2");
  Field f = p.solve();
  ok(f.l2_error("x^2+y^2") < 1e-9, "nonzero polynomial Dirichlet data");
  ok(f.components() == 1 && f.dofs() > 0, "scalar field metadata");
}

static void heat_time_and_boundary() {
  Problem p;
  p.mesh = Mesh::rect(12, 12, 1.0, 1.0);
  p.kind = "heat";
  p.order = 2;
  p.initial("sin(pi*x)*sin(pi*y)").source("0").dirichlet(0, "0");
  Field f = p.evolve(0.01, 0.001);
  double expected = std::exp(-2.0 * M_PI * M_PI * 0.01);
  ok(std::abs(f.value(0.5, 0.5, 0) - expected) < 2e-3, "heat sine mode");
  ok(std::abs(f.time() - 0.01) < 1e-13, "heat final time");

  Problem q;
  q.mesh = Mesh::rect(8, 8, 1.0, 1.0);
  q.kind = "heat";
  q.initial("0").source("1").dirichlet(0, "t");
  Field g = q.evolve(0.02, 0.002);
  ok(std::abs(g.value(0.5, 0.5, 0) - 0.02) < 3e-3, "time-aware heat boundary");
}

static void elasticity_patch() {
  Problem p;
  p.mesh = Mesh::rect(4, 4, 1.0, 1.0);
  p.kind = "elasticity";
  p.order = 1;
  p.lame(2.0, 3.0).displacement(0, "x", "y", "0").force("0", "0", "0");
  Field f = p.solve();
  ok(f.components() == 2, "elasticity component count");
  ok(f.l2_error_vector("x", "y", "0") < 1e-8, "elasticity affine patch test");
}

static void mesh_3d_and_lifetime() {
  const std::string path = "build/fem-roundtrip.mesh";
  Mesh original = Mesh::box(2, 2, 2, 1.0, 1.0, 1.0);
  original.save(path);
  Mesh loaded = Mesh::load(path);
  ok(loaded.dimension() == 3 && loaded.elements() == original.elements(),
     "3D mesh round trip");
  Problem p;
  p.mesh = loaded;
  p.order = 2;
  p.source("-6").dirichlet(0, "x^2+y^2+z^2");
  Field f = p.solve();
  ok(f.components() == 1 &&
         std::abs(f.value(.25, .25, .25) - 3.0 * .25 * .25) < 2e-2,
     "3D Poisson polynomial");
  Field retained;
  {
    Problem q;
    q.mesh = Mesh::rect(3, 3, 1.0, 1.0);
    q.source("0").dirichlet(0, "1");
    retained = q.solve();
  }
  ok(std::abs(retained.value(.5, .5, 0) - 1.0) < 1e-10,
     "field retains mesh ownership");
  std::remove(path.c_str());
}

static void invalid_points() {
  Problem p;
  p.mesh = Mesh::rect(2, 2, 1.0, 1.0);
  p.source("0").dirichlet(0, "0");
  Field f = p.solve();
  ok(throws_bad([&] { (void)f.value(-0.1, 0.5, 0); }),
     "out of domain point rejected");
  ok(throws_bad([] { (void)Mesh::rect(0, 2, 1, 1); }),
     "invalid mesh dimensions rejected");
  ok(throws_bad([] {
       Problem q;
       q.mesh = Mesh::rect(2, 2, 1, 1);
       q.diffusion("0");
       (void)q.solve();
     }),
     "nonpositive coefficient rejected");
}

static void boundary_and_transient_regressions() {
  Problem p;
  p.mesh = Mesh::rect(5, 5, 1, 1);
  p.source("0").dirichlet(4, "0").dirichlet(2, "1");
  ok(p.solve().l2_error("x") < 1e-9,
     "mixed scalar Dirichlet and natural faces");
  Problem e;
  e.mesh = p.mesh;
  e.kind = "elasticity";
  e.lame(0, 1).displacement(4, "0", "0", "0").displacement(2, "1", "0", "0");
  ok(e.solve().l2_error_vector("x", "0", "0") < 1e-9,
     "attribute-specific vector boundary data");
  Problem h;
  h.mesh = p.mesh;
  h.kind = "heat";
  h.order = 2;
  h.initial("0")
      .diffusion("1+t")
      .source("x*x+y*y-4*t*(1+t)")
      .dirichlet(0, "t*(x*x+y*y)");
  ok(h.evolve(.1, .02).l2_error("t*(x*x+y*y)") < 1e-9,
     "time varying coefficient/source/Dirichlet manufactured heat");
  h.initial("0").source("1").dirichlet(0, "t");
  auto tiny = h.evolve(1e-16, 1e-16);
  ok(tiny.time() == 1e-16, "tiny heat end time advanced");
  ok(std::abs(tiny.value(.5, .5, 0) - 1e-16) < 1e-25,
     "tiny heat solution advanced");
  h.dirichlet(0, "2");
  ok(std::abs(h.evolve(0, .1).value(0, .5, 0) - 2) < 1e-14,
     "initial heat boundary projected at t0");
  ok(throws_bad([&] {
       Problem q;
       q.mesh = p.mesh;
       q.diffusion("1-2*x");
       q.solve();
     }),
     "quadrature detects spatially negative diffusion");
  ok(throws_bad([&] {
       Problem q;
       q.mesh = p.mesh;
       q.dirichlet(99, "0").solve();
     }),
     "missing boundary attribute rejected");
  ok(throws_bad([&] {
       Problem q;
       q.mesh = p.mesh;
       q.source("sqrt(-1)").solve();
     }),
     "nonfinite source rejected");
  ok(throws_bad([] { Mesh::rect(1000000, 1000000, 1, 1); }),
     "oversized mesh rejected before allocation");
}

static void higher_order_and_body_force() {
  Problem cubic;
  cubic.mesh = Mesh::rect(3, 3, 1, 1);
  cubic.order = 3;
  cubic.source("-6*x-6*y").dirichlet(0, "x^3+y^3");
  ok(cubic.solve().l2_error("x^3+y^3") < 1e-8, "cubic H1 polynomial patch");
  Problem elastic;
  elastic.mesh = cubic.mesh;
  elastic.kind = "elasticity";
  elastic.order = 2;
  elastic.lame(1, 1).force("-6", "-6", "0").displacement(0, "x*x", "y*y", "0");
  ok(elastic.solve().l2_error_vector("x*x", "y*y", "0") < 1e-8,
     "elasticity nonzero body force patch");
  elastic.mesh = Mesh::box(2, 2, 2, 1, 1, 1);
  elastic.order = 1;
  elastic.force("0", "0", "0").displacement(0, "x", "y", "z");
  ok(elastic.solve().l2_error_vector("x", "y", "z") < 1e-8,
     "3D elasticity affine patch");
}

int main() {
  try {
    convergence();
    scalar_polynomial_and_coefficient();
    heat_time_and_boundary();
    elasticity_patch();
    mesh_3d_and_lifetime();
    invalid_points();
    boundary_and_transient_regressions();
    higher_order_and_body_force();
    std::cout << "fem checks: " << checks << " passed\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "fem failure: " << e.what() << "\n";
    return 1;
  }
}
