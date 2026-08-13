#ifndef G4BRIDGE_RESULTS_HPP
#define G4BRIDGE_RESULTS_HPP

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace g4bridge
{

struct Results
{
  // session state
  bool initialized = false;
  bool has_source = false;
  std::string physics_list = "";
  std::size_t loaded_primaries = 0;

  // most recent run summary
  std::size_t last_events_run = 0;
  double last_total_edep_mev = 0.0;
  double last_dose_gy = 0.0;
  std::vector<double> edep_spectrum_edges_mev;
  std::vector<std::size_t> edep_spectrum_counts;
  std::vector<double> edep_spectrum_edep_mev;
  std::size_t edep_spectrum_underflow = 0;
  std::size_t edep_spectrum_overflow = 0;
  std::vector<std::string> component_names;
  std::vector<double> component_edep_mev;
  std::vector<double> component_mass_kg;
  std::vector<double> component_dose_gy;

  std::string status = "created";

  // clear fields that are produced by a run or source load
  void ResetRuntime();
};

}  // namespace g4bridge

#endif
