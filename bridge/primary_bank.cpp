#include "primary_bank.hpp"

#include <cmath>
#include <sstream>
#include <stdexcept>
#include <utility>

namespace py = pybind11;

namespace g4bridge
{

namespace
{

constexpr py::ssize_t kExpectedCols = 10;
constexpr double kIntegerTolerance = 1e-9;

enum PrimaryCol : py::ssize_t {
  kParticleId = 0,
  kXmm = 1,
  kYmm = 2,
  kZmm = 3,
  kUx = 4,
  kUy = 5,
  kUz = 6,
  kEnergyMeV = 7,
  kWeight = 8,
  kTimeNs = 9
};

void RequireFinite(double value, std::size_t row, const char* name)
{
  if (!std::isfinite(value)) {
    std::ostringstream msg;
    msg << "Invalid primary row " << row << ": '" << name << "' must be finite.";
    throw std::runtime_error(msg.str());
  }
}

int ParseParticleId(double raw_id, std::size_t row)
{
  // convert stored particle id to an integer Geant4 id
  RequireFinite(raw_id, row, "particle_id");
  const auto rounded = static_cast<long long>(std::llround(raw_id));
  if (std::fabs(raw_id - static_cast<double>(rounded)) > kIntegerTolerance) {
    std::ostringstream msg;
    msg << "Invalid primary row " << row << ": particle_id must be an integer value.";
    throw std::runtime_error(msg.str());
  }
  return static_cast<int>(rounded);
}

}  // namespace

void PrimaryBank::Load(
  const py::array_t<double, py::array::c_style | py::array::forcecast>& bank)
{
  // validate primary bank shape
  const py::buffer_info info = bank.request();
  if (info.ndim != 2) {
    throw std::runtime_error("Primary bank must be a 2D NumPy array with shape (N, 10).");
  }
  if (info.shape[1] != kExpectedCols) {
    std::ostringstream msg;
    msg << "Primary bank must have exactly " << kExpectedCols << " columns, got " << info.shape[1]
        << ".";
    throw std::runtime_error(msg.str());
  }

  // copy and validate each primary row
  const auto view = bank.unchecked<2>();
  rows_.clear();
  rows_.reserve(static_cast<std::size_t>(info.shape[0]));

  for (py::ssize_t i = 0; i < info.shape[0]; ++i) {
    Primary p{};
    const auto row = static_cast<std::size_t>(i);

    const double raw_id = view(i, kParticleId);
    p.particle_id = ParseParticleId(raw_id, row);

    p.x_mm = view(i, kXmm);
    p.y_mm = view(i, kYmm);
    p.z_mm = view(i, kZmm);
    p.ux = view(i, kUx);
    p.uy = view(i, kUy);
    p.uz = view(i, kUz);
    p.energy_mev = view(i, kEnergyMeV);
    p.weight = view(i, kWeight);
    p.time_ns = view(i, kTimeNs);

    RequireFinite(p.x_mm, row, "x_mm");
    RequireFinite(p.y_mm, row, "y_mm");
    RequireFinite(p.z_mm, row, "z_mm");
    RequireFinite(p.ux, row, "ux");
    RequireFinite(p.uy, row, "uy");
    RequireFinite(p.uz, row, "uz");
    RequireFinite(p.energy_mev, row, "E_MeV");
    RequireFinite(p.weight, row, "weight");
    RequireFinite(p.time_ns, row, "time_ns");

    if (p.energy_mev < 0.0) {
      std::ostringstream msg;
      msg << "Invalid primary row " << row << ": E_MeV must be >= 0.";
      throw std::runtime_error(msg.str());
    }

    const double norm = std::sqrt(p.ux * p.ux + p.uy * p.uy + p.uz * p.uz);
    if (!(norm > 0.0)) {
      std::ostringstream msg;
      msg << "Invalid primary row " << row << ": direction (ux,uy,uz) norm must be > 0.";
      throw std::runtime_error(msg.str());
    }
    // normalize direction before handing it to Geant4
    p.ux /= norm;
    p.uy /= norm;
    p.uz /= norm;

    rows_.emplace_back(std::move(p));
  }
}

void PrimaryBank::Clear()
{
  rows_.clear();
}

}  // namespace g4bridge
