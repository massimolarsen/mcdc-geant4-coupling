#ifndef G4BRIDGE_SESSION_HPP
#define G4BRIDGE_SESSION_HPP

#include <array>
#include <memory>
#include <string>
#include <vector>

#include <pybind11/numpy.h>

#include "device_geometry.hpp"
#include "results.hpp"

class G4RunManager;
class G4UIsession;

namespace g4bridge
{

class PrimaryBank;
class SourceDistribution;
struct ProgressState;

struct SessionConfig
{
  // detector setup used to construct the Geant4 world
  std::array<double, 3> world_size_mm{100.0, 100.0, 100.0};
  std::array<double, 3> detector_size_mm{10.0, 10.0, 10.0};
  std::string detector_material = "G4_Si";
  std::string envelope_material = "G4_Galactic";
  std::vector<DeviceComponent> device_components;
  std::string physics_list = "QGSP_BIC";
  long random_seed = 1;
  int n_threads = 1;
  double em_production_cut_mm = 0.0;
  bool record_seu_events = false;
  double diagnostic_min_Eion_mev = 0.001;
  std::string diagnostic_dir;
};

class Session
{
  public:
    // create reusable bridge session
    explicit Session(SessionConfig config);
    ~Session();
    Session(const Session&) = delete;
    Session& operator=(const Session&) = delete;
    Session(Session&&) = delete;
    Session& operator=(Session&&) = delete;

    // initialize Geant4 objects
    void initialize();

    // load exact MCDC handoff particles
    void load_primaries(
      const pybind11::array_t<double, pybind11::array::c_style | pybind11::array::forcecast>&
        primaries);

    // append one species' sampled source distribution
    void load_source_distribution(
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

    // drop loaded source-distribution components
    void clear_source_distributions();

    // run the active Geant4 source
    void beam_on();

    // return latest bridge state and run summary
    [[nodiscard]] Results get_results() const;

    // release Geant4 run manager
    void close();

  private:
    // validate geometry and physics settings
    void ValidateGeometryConfig(const SessionConfig& config);
    std::string ResolvePhysicsList(const SessionConfig& config) const;

    bool initialized_ = false;
    bool has_source_ = false;

    SessionConfig config_;
    std::string physics_list_name_;

    std::unique_ptr<G4RunManager> run_manager_;
    std::unique_ptr<G4UIsession> silent_ui_session_;
    std::unique_ptr<PrimaryBank> primary_bank_;
    std::unique_ptr<SourceDistribution> source_distribution_;
    std::shared_ptr<ProgressState> progress_state_;
    Results results_;
};

}  // namespace g4bridge

#endif
