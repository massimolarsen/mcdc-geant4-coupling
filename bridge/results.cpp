#include "results.hpp"

namespace g4bridge
{

void Results::ResetRuntime()
{
  // clear per-run result fields
  last_events_run = 0;
  last_total_edep_mev = 0.0;
  last_dose_gy = 0.0;
  edep_spectrum_edges_mev.clear();
  edep_spectrum_counts.clear();
  edep_spectrum_edep_mev.clear();
  edep_spectrum_underflow = 0;
  edep_spectrum_overflow = 0;
  component_names.clear();
  component_edep_mev.clear();
  component_mass_kg.clear();
  component_dose_gy.clear();
  component_niel_mev.clear();
  component_ionizing_mev.clear();
  component_ionizing_sum_sq_mev2.clear();
  seu_species_names.clear();
  component_species_ionizing_mev.clear();
  component_species_ionizing_sum_sq_mev2.clear();
  component_species_positive_events.clear();
  component_primary_ionizing_mev.clear();
  component_secondary_ionizing_mev.clear();
  component_event_ionizing_edges_mev.clear();
  component_event_ionizing_count.clear();
  component_event_ionizing_sumw.clear();
  component_event_ionizing_sumw2.clear();
  component_primary_species_edep_mev.clear();
  component_primary_species_ionizing_mev.clear();
  component_primary_species_ionizing_sum_sq_mev2.clear();
  component_primary_species_event_ionizing_count.clear();
  component_primary_species_event_ionizing_sumw.clear();
  component_primary_species_event_ionizing_sumw2.clear();
  electronics_cut_materials.clear();
  electronics_cut_energy_mev.clear();
  tracks.Clear();
}

}  // namespace g4bridge
