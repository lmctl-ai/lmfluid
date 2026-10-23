#pragma once
#include <memory>
#include <string>
#include <vector>
namespace lmfluid {
struct Matrix {
  int rows = 0, cols = 0;
  std::vector<double> data;
  double get(int, int) const;
  void save_csv(const std::string &) const;
};
struct MeshData;
struct FieldData;
struct Mesh {
  std::shared_ptr<MeshData> data;
  static Mesh rect(int nx, int ny, double lx, double ly);
  static Mesh box(int nx, int ny, int nz, double lx, double ly, double lz);
  static Mesh load(const std::string &);
  int dimension() const;
  int elements() const;
  void save(const std::string &) const;
};
struct Field {
  std::shared_ptr<FieldData> data;
  double value(double x, double y, double z, int component = 0) const;
  double l2_error(const std::string &) const;
  double l2_error_vector(const std::string &, const std::string &,
                         const std::string &) const;
  Matrix sample(int nx, int ny, int nz = 1) const;
  int dofs() const;
  int components() const;
  double time() const;
};
struct Boundary {
  int attribute = 0;
  std::string x = "0", y = "0", z = "0";
};
struct Problem {
  Mesh mesh;
  std::string kind = "poisson";
  int order = 1;
  std::string coefficient = "1", source_expr = "0", initial_expr = "0";
  std::string force_x = "0", force_y = "0", force_z = "0";
  double lambda = 1, mu = 1, relative_tolerance = 1e-10;
  int max_iterations = 10000;
  std::vector<Boundary> boundaries{{0, "0", "0", "0"}};
  Problem &diffusion(const std::string &);
  Problem &source(const std::string &);
  Problem &initial(const std::string &);
  Problem &dirichlet(int, const std::string &);
  Problem &displacement(int, const std::string &, const std::string &,
                        const std::string &);
  Problem &force(const std::string &, const std::string &, const std::string &);
  Problem &lame(double, double);
  Problem &tol(double, int);
  Field solve() const;
  Field evolve(double t_end, double dt) const;
};
struct LbmData;
struct LBM {
  std::shared_ptr<LbmData> data;
  static LBM channel(int nx, int ny, double viscosity, double acceleration);
  static LBM cavity(int nx, int ny, double viscosity, double lid_speed);
  static LBM periodic(int nx, int ny, double viscosity);
  LBM &init(const std::string &ux, const std::string &uy);
  LBM &run(int steps);
  double velocity(int ix, int iy, int component) const;
  double density(int ix, int iy) const;
  double mass() const;
  double max_speed() const;
  double kinetic_energy() const;
  int steps() const;
  Matrix sample() const;
};
} // namespace lmfluid
