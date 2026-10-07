#ifndef G4BRIDGE_RESULTS_HPP
#define G4BRIDGE_RESULTS_HPP

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "recorded_tracks.hpp"

namespace g4bridge
{

struct Results
{
  // session state
  std::string physics_list = "";
  std::string geant4_version = "";
  std::size_t loaded_primaries = 0;
  double em_production_cut_mm = 0.0;
  double proton_production_cut_mm = 0.0;
  std::vector<std::string> electronics_cut_materials;
  // One row per material; columns are gamma, e-, e+, proton.
  std::vector<double> electronics_cut_energy_mev;

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
  std::vector<double> component_niel_mev;
  std::vector<double> component_ionizing_mev;
  std::vector<double> component_ionizing_sum_sq_mev2;
  std::vector<std::string> seu_species_names;
  std::vector<double> component_species_ionizing_mev;
  std::vector<double> component_species_ionizing_sum_sq_mev2;
  std::vector<std::size_t> component_species_positive_events;
  std::vector<double> component_primary_ionizing_mev;
  std::vector<double> component_secondary_ionizing_mev;
  std::vector<double> component_event_ionizing_edges_mev;
  std::vector<double> component_event_ionizing_count;
  std::vector<double> component_event_ionizing_sumw;
  std::vector<double> component_event_ionizing_sumw2;
  // split by primary species (SeU species order): component x species [x SeU bin]
  std::vector<double> component_primary_species_edep_mev;
  std::vector<double> component_primary_species_ionizing_mev;
  std::vector<double> component_primary_species_ionizing_sum_sq_mev2;
  std::vector<double> component_primary_species_event_ionizing_count;
  std::vector<double> component_primary_species_event_ionizing_sumw;
  std::vector<double> component_primary_species_event_ionizing_sumw2;

  // tracks of the last event, when SessionConfig::record_tracks is set
  RecordedTracks tracks;

  // clear fields that are produced by a run or source load
  void ResetRuntime();
};

}  // namespace g4bridge

#endif
