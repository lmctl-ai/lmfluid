#include "core.hpp"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
namespace lmfluid {
double Matrix::get(int i, int j) const {
  if (i < 0 || j < 0 || i >= rows || j >= cols)
    throw std::runtime_error("Matrix index outside range");
  return data.at(static_cast<size_t>(i) * cols + j);
}
void Matrix::set(int i, int j, double v) {
  if (i < 0 || j < 0 || i >= rows || j >= cols)
    throw std::runtime_error("Matrix index outside range");
  if (!std::isfinite(v))
    throw std::runtime_error("Matrix value must be finite");
  data.at(static_cast<size_t>(i) * cols + j) = v;
}
void Matrix::save_csv(const std::string &path) const {
  if (rows < 0 || cols < 0 || data.size() != static_cast<size_t>(rows) * cols)
    throw std::runtime_error("invalid Matrix shape");
  std::ofstream file(path);
  if (!file)
    throw std::runtime_error("cannot open CSV: " + path);
  file << std::setprecision(17);
  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      double v = get(i, j);
      if (!std::isfinite(v))
        throw std::runtime_error("nonfinite CSV value");
      if (j)
        file << ',';
      file << v;
    }
    file << '\n';
  }
  file.flush();
  if (!file)
    throw std::runtime_error("cannot write CSV: " + path);
}
} // namespace lmfluid
