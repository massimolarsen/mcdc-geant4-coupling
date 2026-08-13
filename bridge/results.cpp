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
  status = "reset";
}

}  // namespace g4bridge
