#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "session.hpp"
#ifdef G4BRIDGE_HAS_VIEWER
#include "event_viewer.hpp"
#endif

namespace py = pybind11;

PYBIND11_MODULE(geant4_bridge, m)
{
  m.doc() = "Geant4 bridge module with reusable session API";

  // expose local device component geometry
  py::class_<g4bridge::DeviceComponent>(m, "DeviceComponent")
    .def(py::init<>())
    .def_readwrite("name", &g4bridge::DeviceComponent::name)
    .def_readwrite("material", &g4bridge::DeviceComponent::material)
    .def_readwrite("parent", &g4bridge::DeviceComponent::parent)
    .def_readwrite("center_mm", &g4bridge::DeviceComponent::center_mm)
    .def_readwrite("size_mm", &g4bridge::DeviceComponent::size_mm)
    .def_readwrite("score", &g4bridge::DeviceComponent::score);

  // expose session configuration
  py::class_<g4bridge::SessionConfig>(m, "SessionConfig")
    .def(py::init<>())
    .def_readwrite("world_size_mm", &g4bridge::SessionConfig::world_size_mm)
    .def_readwrite("detector_size_mm", &g4bridge::SessionConfig::detector_size_mm)
    .def_readwrite("detector_material", &g4bridge::SessionConfig::detector_material)
    .def_readwrite("envelope_material", &g4bridge::SessionConfig::envelope_material)
    .def_readwrite("device_components", &g4bridge::SessionConfig::device_components)
    .def_readwrite("physics_list", &g4bridge::SessionConfig::physics_list)
    .def_readwrite("random_seed", &g4bridge::SessionConfig::random_seed)
    .def_readwrite("n_threads", &g4bridge::SessionConfig::n_threads)
    .def_readwrite("em_production_cut_mm", &g4bridge::SessionConfig::em_production_cut_mm)
    .def_readwrite("record_seu_events", &g4bridge::SessionConfig::record_seu_events)
    .def_readwrite("diagnostic_min_Eion_mev", &g4bridge::SessionConfig::diagnostic_min_Eion_mev)
    .def_readwrite("diagnostic_dir", &g4bridge::SessionConfig::diagnostic_dir)
    .def_readwrite("rng_state_min_Eion_mev", &g4bridge::SessionConfig::rng_state_min_Eion_mev)
    .def_readwrite("replay_mode", &g4bridge::SessionConfig::replay_mode)
    .def_readwrite("event_id_offset", &g4bridge::SessionConfig::event_id_offset)
    .def_readwrite("replay_state", &g4bridge::SessionConfig::replay_state)
    .def_readwrite("record_tracks", &g4bridge::SessionConfig::record_tracks);

  // expose one event's tracks and step points
  py::class_<g4bridge::RecordedTracks>(m, "RecordedTracks")
    .def(py::init<>())
    .def_readwrite("track_id", &g4bridge::RecordedTracks::track_id)
    .def_readwrite("parent_id", &g4bridge::RecordedTracks::parent_id)
    .def_readwrite("pdg", &g4bridge::RecordedTracks::pdg)
    .def_readwrite("particle", &g4bridge::RecordedTracks::particle)
    .def_readwrite("charge", &g4bridge::RecordedTracks::charge)
    .def_readwrite("creator_process", &g4bridge::RecordedTracks::creator_process)
    .def_readwrite("track_point_start", &g4bridge::RecordedTracks::track_point_start)
    .def_readwrite("x_mm", &g4bridge::RecordedTracks::x_mm)
    .def_readwrite("y_mm", &g4bridge::RecordedTracks::y_mm)
    .def_readwrite("z_mm", &g4bridge::RecordedTracks::z_mm)
    .def_readwrite("t_ns", &g4bridge::RecordedTracks::t_ns)
    .def_readwrite("ke_mev", &g4bridge::RecordedTracks::ke_mev)
    .def_readwrite("edep_mev", &g4bridge::RecordedTracks::edep_mev)
    .def_readwrite("process", &g4bridge::RecordedTracks::process)
    .def_readwrite("volume", &g4bridge::RecordedTracks::volume);

  // expose run results
  py::class_<g4bridge::Results>(m, "RunResult")
    .def_readonly("physics_list", &g4bridge::Results::physics_list)
    .def_readonly("loaded_primaries", &g4bridge::Results::loaded_primaries)
    .def_readonly("last_events_run", &g4bridge::Results::last_events_run)
    .def_readonly("last_total_edep_mev", &g4bridge::Results::last_total_edep_mev)
    .def_readonly("last_dose_gy", &g4bridge::Results::last_dose_gy)
    .def_readonly("edep_spectrum_edges_mev", &g4bridge::Results::edep_spectrum_edges_mev)
    .def_readonly("edep_spectrum_counts", &g4bridge::Results::edep_spectrum_counts)
    .def_readonly("edep_spectrum_edep_mev", &g4bridge::Results::edep_spectrum_edep_mev)
    .def_readonly("edep_spectrum_underflow", &g4bridge::Results::edep_spectrum_underflow)
    .def_readonly("edep_spectrum_overflow", &g4bridge::Results::edep_spectrum_overflow)
    .def_readonly("component_names", &g4bridge::Results::component_names)
    .def_readonly("component_edep_mev", &g4bridge::Results::component_edep_mev)
    .def_readonly("component_mass_kg", &g4bridge::Results::component_mass_kg)
    .def_readonly("component_dose_gy", &g4bridge::Results::component_dose_gy)
    .def_readonly("geant4_version", &g4bridge::Results::geant4_version)
    .def_readonly("em_production_cut_mm", &g4bridge::Results::em_production_cut_mm)
    .def_readonly("proton_production_cut_mm", &g4bridge::Results::proton_production_cut_mm)
    .def_readonly("electronics_cut_materials", &g4bridge::Results::electronics_cut_materials)
    .def_readonly("electronics_cut_energy_mev", &g4bridge::Results::electronics_cut_energy_mev)
    .def_readonly("component_niel_mev", &g4bridge::Results::component_niel_mev)
    .def_readonly("component_ionizing_mev", &g4bridge::Results::component_ionizing_mev)
    .def_readonly("component_ionizing_sum_sq_mev2", &g4bridge::Results::component_ionizing_sum_sq_mev2)
    .def_readonly("seu_species_names", &g4bridge::Results::seu_species_names)
    .def_readonly("component_species_ionizing_mev", &g4bridge::Results::component_species_ionizing_mev)
    .def_readonly("component_species_ionizing_sum_sq_mev2", &g4bridge::Results::component_species_ionizing_sum_sq_mev2)
    .def_readonly("component_species_positive_events", &g4bridge::Results::component_species_positive_events)
    .def_readonly("component_primary_ionizing_mev", &g4bridge::Results::component_primary_ionizing_mev)
    .def_readonly("component_secondary_ionizing_mev", &g4bridge::Results::component_secondary_ionizing_mev)
    .def_readonly("component_event_ionizing_edges_mev", &g4bridge::Results::component_event_ionizing_edges_mev)
    .def_readonly("component_event_ionizing_count", &g4bridge::Results::component_event_ionizing_count)
    .def_readonly("component_event_ionizing_sumw", &g4bridge::Results::component_event_ionizing_sumw)
    .def_readonly("component_event_ionizing_sumw2", &g4bridge::Results::component_event_ionizing_sumw2)
    .def_readonly(
      "component_primary_species_edep_mev",
      &g4bridge::Results::component_primary_species_edep_mev)
    .def_readonly(
      "component_primary_species_ionizing_mev",
      &g4bridge::Results::component_primary_species_ionizing_mev)
    .def_readonly(
      "component_primary_species_ionizing_sum_sq_mev2",
      &g4bridge::Results::component_primary_species_ionizing_sum_sq_mev2)
    .def_readonly(
      "component_primary_species_event_ionizing_count",
      &g4bridge::Results::component_primary_species_event_ionizing_count)
    .def_readonly(
      "component_primary_species_event_ionizing_sumw",
      &g4bridge::Results::component_primary_species_event_ionizing_sumw)
    .def_readonly(
      "component_primary_species_event_ionizing_sumw2",
      &g4bridge::Results::component_primary_species_event_ionizing_sumw2)
    .def_readonly("tracks", &g4bridge::Results::tracks);

  // expose reusable Geant4 session
  py::class_<g4bridge::Session>(m, "Session")
    .def(py::init<g4bridge::SessionConfig>(), py::arg("config"))
    .def("initialize", &g4bridge::Session::initialize)
    .def("load_primaries", &g4bridge::Session::load_primaries)
    .def(
      "load_source_distribution",
      &g4bridge::Session::load_source_distribution,
      py::arg("box_bounds_mm"),
      py::arg("mu_edges"),
      py::arg("azi_edges"),
      py::arg("energy_edges_mev"),
      py::arg("weights"),
      py::arg("n_u"),
      py::arg("n_v"),
      py::arg("n_events"),
      py::arg("particle_id") = 2112)
    .def("clear_source_distributions", &g4bridge::Session::clear_source_distributions)
    .def(
      "beam_on",
      &g4bridge::Session::beam_on,
      py::arg("n_events") = -1,
      py::call_guard<py::gil_scoped_release>())
    .def("get_results", &g4bridge::Session::get_results)
    .def("close", &g4bridge::Session::close);

  m.attr("has_viewer") = py::bool_(
#ifdef G4BRIDGE_HAS_VIEWER
    true
#else
    false
#endif
  );
#ifdef G4BRIDGE_HAS_VIEWER
  m.def(
    "view_event",
    &g4bridge::ViewEvent,
    py::arg("config"),
    py::arg("tracks"),
    py::arg("commands"),
    "Open Geant4's Qt viewer on recorded tracks inside the config geometry.");
#endif
}
