#include "core.hpp"
#include "dslvm.hpp"
#include "expr.hpp"
#include "memory_budget.hpp"
#include "notices_embedded.hpp"
#include "schema_embedded.hpp"
#include <cmath>
#include <cstdio>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#ifdef _WIN32
#include <fcntl.h>
#include <io.h>
#endif
namespace {
using dslvm::Args;
using dslvm::Value;
using lmfluid::Field;
using lmfluid::LBM;
using lmfluid::Matrix;
using lmfluid::Mesh;
using lmfluid::Problem;
struct Fluid {};
struct Host {
  std::vector<std::unique_ptr<Mesh>> meshes;
  std::vector<std::unique_ptr<Problem>> problems;
  std::vector<std::unique_ptr<Field>> fields;
  std::vector<std::unique_ptr<LBM>> lbms;
  std::vector<std::unique_ptr<Matrix>> matrices;
  template <class T>
  Value hold(std::vector<std::unique_ptr<T>> &arena, const char *type,
             T value) {
    if (arena.size() >= 100000)
      throw dslvm::Error(std::string(type) + " object limit exceeded");
    arena.push_back(std::make_unique<T>(std::move(value)));
    return Value::object(type, arena.back().get());
  }
  Value mesh(Mesh m) { return hold(meshes, "Mesh", std::move(m)); }
  Value field(Field f) { return hold(fields, "Field", std::move(f)); }
  Value lbm(LBM f) { return hold(lbms, "LBM", std::move(f)); }
  Value matrix(Matrix m) { return hold(matrices, "Matrix", std::move(m)); }
  Value problem(const std::string &kind, const Mesh &mesh, int order) {
    if (order < 1 || order > 3)
      throw dslvm::Error("FEM order must be in [1,3]");
    Problem p;
    p.kind = kind;
    p.mesh = mesh;
    p.order = order;
    return hold(problems, "Problem", std::move(p));
  }
};
template <class F> Value translate(F &&f) {
  try {
    return f();
  } catch (const dslvm::Error &) {
    throw;
  } catch (const std::exception &e) {
    // MFEM errors append "\n ... in function: ...\n ... in file: <build path>:<line>"; keep the message
    std::string msg = e.what();
    const auto trailer = msg.find("\n ... in function");
    if (trailer != std::string::npos)
      msg.erase(trailer);
    throw dslvm::Error(msg);
  }
}

void check(dslvm::Registry &r, const char *type, const char *name,
           std::vector<std::string> ps, std::string ret, dslvm::Method fn) {
  const std::string who = std::string(type) + "." + name;
  const bool real = ret == "float";
  r.method(type, name, ps, ret, [who, ps, fn, real](void *self, Args &a) {
    for (size_t i = 0; i < ps.size(); ++i) {
      const auto &p = a[i];
      const bool ok = ps[i] == "float"    ? p.kind == Value::FLOAT
                      : ps[i] == "int"    ? p.kind == Value::INT
                      : ps[i] == "string" ? p.kind == Value::STR
                      : ps[i] == "bool"
                          ? p.kind == Value::BOOL
                          : p.kind == Value::OBJ && p.obj && p.type == ps[i];
      if (!ok)
        throw dslvm::Error(who + " argument " + std::to_string(i + 1) +
                           " must be " + ps[i]);
    }
    Value out = translate([&] { return fn(self, a); });
    // never hand NaN/inf to a script with ok: true
    if (real && out.kind == Value::FLOAT && !std::isfinite(out.f))
      throw dslvm::Error(who + ": the result is not finite");
    return out;
  });
}

int integer(const Value &v, int lo, int hi, const char *what) {
  if (v.i < lo || v.i > hi)
    throw dslvm::Error(std::string(what) + " outside allowed range");
  return static_cast<int>(v.i);
}
Mesh &mesh(void *p) { return *static_cast<Mesh *>(p); }
Problem &problem(void *p) { return *static_cast<Problem *>(p); }
Field &field(void *p) { return *static_cast<Field *>(p); }
LBM &lbm(void *p) { return *static_cast<LBM *>(p); }
Matrix &matrix(void *p) { return *static_cast<Matrix *>(p); }
std::string fmt(double x) {
  char b[48];
  std::snprintf(b, sizeof b, "%.17g", x);
  return b;
}
void register_host(dslvm::Registry &r, Host &h) {
  r.context("fluid", "Fluid");
  check(r, "Fluid", "rect", {"int", "int", "float", "float"}, "Mesh",
        [&h](void *, Args &a) {
          return h.mesh(Mesh::rect(integer(a[0], 1, 4096, "mesh nx"),
                                   integer(a[1], 1, 4096, "mesh ny"), a[2].f,
                                   a[3].f));
        });
  check(r, "Fluid", "box", {"int", "int", "int", "float", "float", "float"},
        "Mesh", [&h](void *, Args &a) {
          return h.mesh(Mesh::box(integer(a[0], 1, 512, "mesh nx"),
                                  integer(a[1], 1, 512, "mesh ny"),
                                  integer(a[2], 1, 512, "mesh nz"), a[3].f,
                                  a[4].f, a[5].f));
        });
  check(r, "Fluid", "load_mesh", {"string"}, "Mesh",
        [&h](void *, Args &a) { return h.mesh(Mesh::load(a[0].s)); });
  check(r, "Fluid", "poisson", {"Mesh", "int"}, "Problem",
        [&h](void *, Args &a) {
          return h.problem("poisson", mesh(a[0].obj),
                           integer(a[1], 1, 3, "FEM order"));
        });
  check(r, "Fluid", "heat", {"Mesh", "int"}, "Problem", [&h](void *, Args &a) {
    return h.problem("heat", mesh(a[0].obj), integer(a[1], 1, 3, "FEM order"));
  });
  check(r, "Fluid", "elasticity", {"Mesh", "int"}, "Problem",
        [&h](void *, Args &a) {
          return h.problem("elasticity", mesh(a[0].obj),
                           integer(a[1], 1, 3, "FEM order"));
        });
  check(r, "Fluid", "channel", {"int", "int", "float", "float"}, "LBM",
        [&h](void *, Args &a) {
          return h.lbm(LBM::channel(integer(a[0], 3, 2048, "LBM nx"),
                                    integer(a[1], 3, 2048, "LBM ny"), a[2].f,
                                    a[3].f));
        });
  check(r, "Fluid", "cavity", {"int", "int", "float", "float"}, "LBM",
        [&h](void *, Args &a) {
          return h.lbm(LBM::cavity(integer(a[0], 3, 2048, "LBM nx"),
                                   integer(a[1], 3, 2048, "LBM ny"), a[2].f,
                                   a[3].f));
        });
  check(r, "Fluid", "periodic", {"int", "int", "float"}, "LBM",
        [&h](void *, Args &a) {
          return h.lbm(LBM::periodic(integer(a[0], 3, 2048, "LBM nx"),
                                     integer(a[1], 3, 2048, "LBM ny"), a[2].f));
        });
  check(r, "Fluid", "assert_close", {"float", "float", "float"}, "void",
        [](void *, Args &a) {
          if (!std::isfinite(a[0].f) || !std::isfinite(a[1].f) ||
              !std::isfinite(a[2].f) || a[2].f < 0 ||
              std::abs(a[0].f - a[1].f) > a[2].f)
            throw dslvm::Error("assert_close failed: " + fmt(a[0].f) + " vs " +
                               fmt(a[1].f));
          return Value::none();
        });
  check(r, "Fluid", "assert_true", {"bool"}, "void", [](void *, Args &a) {
    if (!a[0].b)
      throw dslvm::Error("assert_true failed");
    return Value::none();
  });
  check(r, "Fluid", "eval", {"string"}, "float", [](void *, Args &a) {
    lmfluid::Expr e(a[0].s, {});
    double x = e.eval(nullptr);
    if (!std::isfinite(x))
      throw dslvm::Error("expression result must be finite");
    return Value::real(x);
  });
  check(r, "Fluid", "print", {"string"}, "void", [](void *, Args &a) {
    std::cout << a[0].s << '\n';
    return Value::none();
  });
  check(r, "Mesh", "dimension", {}, "int",
        [](void *p, Args &) { return Value::integer(mesh(p).dimension()); });
  check(r, "Mesh", "elements", {}, "int",
        [](void *p, Args &) { return Value::integer(mesh(p).elements()); });
  check(r, "Mesh", "save", {"string"}, "void", [](void *p, Args &a) {
    mesh(p).save(a[0].s);
    return Value::none();
  });
  check(r, "Problem", "diffusion", {"string"}, "Problem", [](void *p, Args &a) {
    problem(p).diffusion(a[0].s);
    return Value::object("Problem", p);
  });
  check(r, "Problem", "source", {"string"}, "Problem", [](void *p, Args &a) {
    problem(p).source(a[0].s);
    return Value::object("Problem", p);
  });
  check(r, "Problem", "initial", {"string"}, "Problem", [](void *p, Args &a) {
    problem(p).initial(a[0].s);
    return Value::object("Problem", p);
  });
  check(r, "Problem", "dirichlet", {"int", "string"}, "Problem",
        [](void *p, Args &a) {
          problem(p).dirichlet(integer(a[0], 0, 100000, "boundary attribute"),
                               a[1].s);
          return Value::object("Problem", p);
        });
  check(r, "Problem", "displacement", {"int", "string", "string", "string"},
        "Problem", [](void *p, Args &a) {
          problem(p).displacement(
              integer(a[0], 0, 100000, "boundary attribute"), a[1].s, a[2].s,
              a[3].s);
          return Value::object("Problem", p);
        });
  check(r, "Problem", "force", {"string", "string", "string"}, "Problem",
        [](void *p, Args &a) {
          problem(p).force(a[0].s, a[1].s, a[2].s);
          return Value::object("Problem", p);
        });
  check(r, "Problem", "lame", {"float", "float"}, "Problem",
        [](void *p, Args &a) {
          problem(p).lame(a[0].f, a[1].f);
          return Value::object("Problem", p);
        });
  check(r, "Problem", "tol", {"float", "int"}, "Problem", [](void *p, Args &a) {
    problem(p).tol(a[0].f, integer(a[1], 1, 100000, "solver iterations"));
    return Value::object("Problem", p);
  });
  check(r, "Problem", "solve", {}, "Field",
        [&h](void *p, Args &) { return h.field(problem(p).solve()); });
  check(r, "Problem", "evolve", {"float", "float"}, "Field",
        [&h](void *p, Args &a) {
          return h.field(problem(p).evolve(a[0].f, a[1].f));
        });
  check(r, "Field", "value", {"float", "float"}, "float", [](void *p, Args &a) {
    return Value::real(field(p).value(a[0].f, a[1].f, 0, 0));
  });
  check(r, "Field", "value", {"float", "float", "float", "int"}, "float",
        [](void *p, Args &a) {
          return Value::real(field(p).value(
              a[0].f, a[1].f, a[2].f, integer(a[3], 0, 2, "field component")));
        });
  check(r, "Field", "l2_error", {"string"}, "float", [](void *p, Args &a) {
    return Value::real(field(p).l2_error(a[0].s));
  });
  check(r, "Field", "l2_error_vector", {"string", "string", "string"}, "float",
        [](void *p, Args &a) {
          return Value::real(field(p).l2_error_vector(a[0].s, a[1].s, a[2].s));
        });
  check(r, "Field", "sample", {"int", "int"}, "Matrix", [&h](void *p, Args &a) {
    return h.matrix(field(p).sample(integer(a[0], 2, 4096, "sample nx"),
                                    integer(a[1], 2, 4096, "sample ny"), 1));
  });
  check(r, "Field", "sample", {"int", "int", "int"}, "Matrix",
        [&h](void *p, Args &a) {
          return h.matrix(field(p).sample(integer(a[0], 2, 4096, "sample nx"),
                                          integer(a[1], 2, 4096, "sample ny"),
                                          integer(a[2], 1, 4096, "sample nz")));
        });
  check(r, "Field", "dofs", {}, "int",
        [](void *p, Args &) { return Value::integer(field(p).dofs()); });
  check(r, "Field", "components", {}, "int",
        [](void *p, Args &) { return Value::integer(field(p).components()); });
  check(r, "Field", "time", {}, "float",
        [](void *p, Args &) { return Value::real(field(p).time()); });
  check(r, "LBM", "init", {"string", "string"}, "LBM", [](void *p, Args &a) {
    lbm(p).init(a[0].s, a[1].s);
    return Value::object("LBM", p);
  });
  check(r, "LBM", "run", {"int"}, "LBM", [](void *p, Args &a) {
    lbm(p).run(integer(a[0], 0, 10000000, "LBM steps"));
    return Value::object("LBM", p);
  });
  check(r, "LBM", "velocity", {"int", "int", "int"}, "float",
        [](void *p, Args &a) {
          return Value::real(
              lbm(p).velocity(integer(a[0], 0, 2047, "LBM ix"),
                              integer(a[1], 0, 2047, "LBM iy"),
                              integer(a[2], 0, 1, "LBM component")));
        });
  check(r, "LBM", "density", {"int", "int"}, "float", [](void *p, Args &a) {
    return Value::real(lbm(p).density(integer(a[0], 0, 2047, "LBM ix"),
                                      integer(a[1], 0, 2047, "LBM iy")));
  });
  check(r, "LBM", "mass", {}, "float",
        [](void *p, Args &) { return Value::real(lbm(p).mass()); });
  check(r, "LBM", "max_speed", {}, "float",
        [](void *p, Args &) { return Value::real(lbm(p).max_speed()); });
  check(r, "LBM", "kinetic_energy", {}, "float",
        [](void *p, Args &) { return Value::real(lbm(p).kinetic_energy()); });
  check(r, "LBM", "steps", {}, "int",
        [](void *p, Args &) { return Value::integer(lbm(p).steps()); });
  check(r, "LBM", "sample", {}, "Matrix",
        [&h](void *p, Args &) { return h.matrix(lbm(p).sample()); });
  check(r, "Matrix", "rows", {}, "int",
        [](void *p, Args &) { return Value::integer(matrix(p).rows); });
  check(r, "Matrix", "cols", {}, "int",
        [](void *p, Args &) { return Value::integer(matrix(p).cols); });
  check(r, "Matrix", "get", {"int", "int"}, "float", [](void *p, Args &a) {
    return Value::real(matrix(p).get(integer(a[0], 0, 100000000, "matrix row"),
                                     integer(a[1], 0, 1000, "matrix column")));
  });
  check(r, "Matrix", "save_csv", {"string"}, "void", [](void *p, Args &a) {
    matrix(p).save_csv(a[0].s);
    return Value::none();
  });
}
const char *usage =
    "usage: lmfluid --version | --notices | --dsl-schema "
    "[--json] | --dsl-signatures | [--dry-run] --dsl FILE [FILE ...]\n"
    "environment: LMFLUID_MAX_MEMORY_MB  total memory budget in MiB (default "
    "4096); a program that needs more\n"
    "             stops with a run error (exit 3) instead of allocating\n";
} // namespace

int main(int argc, char **argv) {
#ifdef _WIN32
  _setmode(_fileno(stdout), _O_BINARY);
  _setmode(_fileno(stderr), _O_BINARY);
#endif
  lmfluid::route_mfem_output();
  Host h;
  dslvm::Registry reg;
  register_host(reg, h);
  if (argc == 2 && std::string(argv[1]) == "--version") {
    std::cout << "lmfluid " LMFLUID_VERSION "\n";
    return 0;
  }
  if (argc == 2 && std::string(argv[1]) == "--notices") {
    std::cout.write(reinterpret_cast<const char *>(lmfluid::embedded_notices),
                    (std::streamsize)lmfluid::embedded_notices_size);
    return 0;
  }
  if (argc >= 2 && std::string(argv[1]) == "--dsl-schema") {
    if (argc == 2) {
      for (const auto &s : reg.signatures())
        std::cout << s << '\n';
      return 0;
    }
    if (argc == 3 && std::string(argv[2]) == "--json") {
      std::cout.write(reinterpret_cast<const char *>(lmfluid::embedded_schema),
                      (std::streamsize)lmfluid::embedded_schema_size);
      return 0;
    }
    std::cerr << usage;
    return 2;
  }
  if (argc == 2 && std::string(argv[1]) == "--dsl-signatures") {
    for (const auto &s : reg.signatures())
      std::cout << s << '\n';
    return 0;
  }
  bool dry = false;
  int dsl_arg = 1;
  if (argc > 1 && std::string(argv[dsl_arg]) == "--dry-run") {
    dry = true;
    ++dsl_arg;
  }
  if (argc <= dsl_arg || std::string(argv[dsl_arg]) != "--dsl") {
    std::cerr << usage;
    return 2;
  }
  std::vector<std::string> files;
  for (int i = dsl_arg + 1; i < argc; i++) {
    if (std::string(argv[i]) == "--dry-run") {
      if (dry) {
        std::cerr << usage;
        return 2;
      }
      dry = true;
    } else if (std::string(argv[i]).empty() ||
               std::string(argv[i]).rfind("--", 0) == 0) {
      std::cerr << usage;
      return 2;
    } else
      files.push_back(argv[i]);
  }
  if (files.empty()) {
    std::cerr << usage;
    return 2;
  }
  if (!membudget::configure("lmfluid", "LMFLUID_MAX_MEMORY_MB")) {
    std::cerr << usage;
    return 2;
  }
  const std::string schema(
      reinterpret_cast<const char *>(lmfluid::embedded_schema),
      lmfluid::embedded_schema_size);
  dslvm::Context ctx;
  Fluid fluid;
  ctx.set("fluid", &fluid);
  std::vector<dslvm::Program> programs;
  for (const auto &file : files)
    try {
      programs.push_back(
          dslvm::Program::load(dslvm::read_file(file), schema, reg));
    } catch (const dslvm::Error &e) {
      return dslvm::report_error("load", e);
    } catch (const std::exception &e) {
      return dslvm::report_error("load", dslvm::Error(e.what()));
    }
  if (dry) {
    dslvm::report_ok();
    return 0;
  }
  for (const auto &p : programs) {
    try {
      p.run(ctx, reg);
    } catch (const dslvm::Error &e) {
      return dslvm::report_error("run", e);
    } catch (const std::exception &e) { // memory budget inside the VM
      return dslvm::report_error("run", dslvm::Error(e.what()));
    }
    dslvm::report_ok();
  }
  return 0;
}
