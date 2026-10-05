#ifndef G4BRIDGE_SOURCE_DISTRIBUTION_HPP
#define G4BRIDGE_SOURCE_DISTRIBUTION_HPP

#include <array>
#include <cstddef>
#include <vector>

#include <pybind11/numpy.h>

#include "primary_bank.hpp"

namespace g4bridge
{

// binned source distribution for one particle species
class SourceComponent
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
      std::size_t n_events,
      int particle_id);

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
    int particle_id_ = 2112;
    std::size_t n_events_ = 0;
    std::size_t n_mu_ = 0;
    std::size_t n_azi_ = 0;
    std::size_t n_energy_ = 0;
    std::size_t n_u_ = 0;
    std::size_t n_v_ = 0;
    double total_weight_ = 0.0;
};

// per-species source components, each run for its own fixed event count
class SourceDistribution
{
  public:
    // append one species component; events follow the load order
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
      std::size_t n_events,
      int particle_id = 2112);

    // clear all loaded components
    void Clear();
    bool Loaded() const { return !components_.empty(); }
    std::size_t Size() const { return n_events_total_; }

    // sample the primary for one event from the component owning that event
    Primary Sample(std::size_t event_id) const;

  private:
    std::vector<SourceComponent> components_;
    std::vector<std::size_t> event_ends_;
    std::size_t n_events_total_ = 0;
};

}  // namespace g4bridge

#endif
