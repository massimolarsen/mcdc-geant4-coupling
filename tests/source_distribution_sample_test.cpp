#include "source_distribution.hpp"

#include <cmath>
#include <initializer_list>
#include <iostream>
#include <stdexcept>

#include <pybind11/embed.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

namespace
{

void Require(bool condition, const char* message)
{
  if (!condition) {
    throw std::runtime_error(message);
  }
}

void RequireRange(double value, double lower, double upper, const char* message)
{
  Require(value >= lower && value <= upper, message);
}

py::array_t<double> Array1D(std::initializer_list<double> values)
{
  py::array_t<double> array(values.size());
  auto view = array.mutable_unchecked<1>();
  std::size_t i = 0;
  for (double value : values) {
    view(static_cast<py::ssize_t>(i)) = value;
    ++i;
  }
  return array;
}

}  // namespace

int main()
{
  try {
    py::scoped_interpreter python{};

    py::array_t<double> box_bounds_mm({3, 2});
    auto box = box_bounds_mm.mutable_unchecked<2>();
    box(0, 0) = -10.0;
    box(0, 1) = 10.0;
    box(1, 0) = -20.0;
    box(1, 1) = 20.0;
    box(2, 0) = -30.0;
    box(2, 1) = 30.0;

    for (std::size_t face = 0; face < 6; ++face) {
      py::array_t<double> weights(6 * 2 * 3);
      auto weight = weights.mutable_unchecked<1>();
      for (py::ssize_t i = 0; i < weight.shape(0); ++i) {
        weight(i) = 0.0;
      }

      // C-order shape is (mu, azi, energy, face, u, v); choose high-u, high-v.
      const std::size_t i_u = 1;
      const std::size_t i_v = 2;
      const std::size_t flat = ((face * 2 + i_u) * 3 + i_v);
      weight(static_cast<py::ssize_t>(flat)) = 1.0;

      g4bridge::SourceDistribution distribution;
      distribution.Load(
        box_bounds_mm,
        Array1D({-0.001, 0.001}),
        Array1D({-0.001, 0.001}),
        Array1D({13.99, 14.01}),
        weights,
        2,
        3,
        1);

      const g4bridge::Primary p = distribution.Sample();

      if (face == 0) {
        Require(p.x_mm == -10.0, "xmin sample should be pinned on x_min.");
        RequireRange(p.y_mm, 0.0, 20.0, "xmin sample should use y as u.");
        RequireRange(p.z_mm, 10.0, 30.0, "xmin sample should use z as v.");
      } else if (face == 1) {
        Require(p.x_mm == 10.0, "xmax sample should be pinned on x_max.");
        RequireRange(p.y_mm, 0.0, 20.0, "xmax sample should use y as u.");
        RequireRange(p.z_mm, 10.0, 30.0, "xmax sample should use z as v.");
      } else if (face == 2) {
        Require(p.y_mm == -20.0, "ymin sample should be pinned on y_min.");
        RequireRange(p.x_mm, 0.0, 10.0, "ymin sample should use x as u.");
        RequireRange(p.z_mm, 10.0, 30.0, "ymin sample should use z as v.");
      } else if (face == 3) {
        Require(p.y_mm == 20.0, "ymax sample should be pinned on y_max.");
        RequireRange(p.x_mm, 0.0, 10.0, "ymax sample should use x as u.");
        RequireRange(p.z_mm, 10.0, 30.0, "ymax sample should use z as v.");
      } else if (face == 4) {
        Require(p.z_mm == -30.0, "zmin sample should be pinned on z_min.");
        RequireRange(p.x_mm, 0.0, 10.0, "zmin sample should use x as u.");
        RequireRange(p.y_mm, 20.0 / 3.0, 20.0, "zmin sample should use y as v.");
      } else {
        Require(p.z_mm == 30.0, "zmax sample should be pinned on z_max.");
        RequireRange(p.x_mm, 0.0, 10.0, "zmax sample should use x as u.");
        RequireRange(p.y_mm, 20.0 / 3.0, 20.0, "zmax sample should use y as v.");
      }

      RequireRange(p.energy_mev, 13.99, 14.01, "sample energy is outside bin.");
      Require(std::abs(p.uz) <= 0.001, "sample direction should use selected mu bin.");
    }

    return 0;
  }
  catch (const std::exception& exc) {
    std::cerr << exc.what() << '\n';
    return 1;
  }
}
