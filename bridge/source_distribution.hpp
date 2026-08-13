#ifndef G4BRIDGE_SOURCE_DISTRIBUTION_HPP
#define G4BRIDGE_SOURCE_DISTRIBUTION_HPP

#include <array>
#include <cstddef>
#include <vector>

#include <pybind11/numpy.h>

#include "primary_bank.hpp"

namespace g4bridge
{

class SourceDistribution
{
  public:
    // load binned source distribution from Python
    void Load(
      const pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast>&
        box_bounds_mm,
      const pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast>&
        mu_edges,
      const pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast>&
        azi_edges,
      const pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast>&
        energy_edges_mev,
      const pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast>&
        weights,
      std::size_t n_u,
      std::size_t n_v,
      std::size_t n_events);

    // clear loaded source distribution
    void Clear();
    bool Loaded() const { return loaded_; }
    std::size_t Size() const { return n_events_; }

    // sample one primary particle from the distribution
    Primary Sample() const;

  private:
    // sampled source geometry and bin data
    std::array<std::array<double, 2>, 3> box_bounds_mm_{};
    std::vector<double> mu_edges_;
    std::vector<double> azi_edges_;
    std::vector<double> energy_edges_mev_;
    std::vector<double> cdf_;

    // loaded source dimensions and normalization
    std::size_t n_events_ = 0;
    std::size_t n_mu_ = 0;
    std::size_t n_azi_ = 0;
    std::size_t n_energy_ = 0;
    std::size_t n_u_ = 0;
    std::size_t n_v_ = 0;
    double total_weight_ = 0.0;
    bool loaded_ = false;
};

}  // namespace g4bridge

#endif
