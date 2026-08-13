#ifndef G4BRIDGE_PRIMARY_BANK_HPP
#define G4BRIDGE_PRIMARY_BANK_HPP

#include <cstddef>
#include <vector>

#include <pybind11/numpy.h>

namespace g4bridge
{

struct Primary
{
  // primary particle data passed from MCDC to Geant4
  int particle_id = 0;
  double x_mm = 0.0;
  double y_mm = 0.0;
  double z_mm = 0.0;
  double ux = 0.0;
  double uy = 0.0;
  double uz = 1.0;
  double energy_mev = 0.0;
  double weight = 1.0;
  double time_ns = 0.0;
};

class PrimaryBank
{
  public:
    // load exact primary rows from Python
    void Load(
      const pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast>& bank);
    void Clear();
    std::size_t Size() const { return rows_.size(); }
    bool Empty() const { return rows_.empty(); }
    const Primary& At(std::size_t idx) const { return rows_.at(idx); }

  private:
    // validated primary rows
    std::vector<Primary> rows_;
};

}  // namespace g4bridge

#endif
