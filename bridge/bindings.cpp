#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "session.hpp"

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
    .def_readwrite("n_threads", &g4bridge::SessionConfig::n_threads);

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
    .def_readonly("component_dose_gy", &g4bridge::Results::component_dose_gy);

  // expose reusable Geant4 session
  py::class_<g4bridge::Session>(m, "Session")
    .def(py::init<g4bridge::SessionConfig>(), py::arg("config"))
    .def("initialize", &g4bridge::Session::initialize)
    .def("load_primaries", &g4bridge::Session::load_primaries)
    .def("load_source_distribution", &g4bridge::Session::load_source_distribution)
    .def("beam_on", &g4bridge::Session::beam_on, py::call_guard<py::gil_scoped_release>())
    .def("get_results", &g4bridge::Session::get_results)
    .def("close", &g4bridge::Session::close);
}
