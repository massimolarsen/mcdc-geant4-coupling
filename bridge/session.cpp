#include "session.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include "bank_primary_generator_action.hpp"
#include "detector_construction.hpp"
#include "primary_bank.hpp"
#include "source_distribution.hpp"

#include "G4Accumulable.hh"
#include "G4AccumulableManager.hh"
#include "G4AccVector.hh"
#include "G4EmParameters.hh"
#include "G4Event.hh"
#include "G4HadronicParameters.hh"
#include "G4LogicalVolume.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#endif
#include "G4PhysListFactory.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4Step.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "G4UIsession.hh"
#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4VUserActionInitialization.hh"
#include "CLHEP/Random/Random.h"

namespace py = pybind11;

namespace g4bridge
{

namespace
{

constexpr std::size_t kEdepSpectrumBins = 20;
constexpr G4double kEdepSpectrumMinMeV = 1.0e-12;
constexpr G4double kEdepSpectrumMaxMeV = 20.0;
constexpr double kGeometryToleranceMm = 1.0e-9;

std::size_t CountScoringComponents(const std::vector<DeviceComponent>& components)
{
  return static_cast<std::size_t>(
    std::count_if(components.begin(), components.end(), [](const auto& component) {
      return component.score;
    }));
}

class BridgeRunAction : public G4UserRunAction
{
  public:
    BridgeRunAction(Results& results, std::size_t n_components, bool record_results)
        : results_(results),
          record_results_(record_results),
          component_edep_(n_components),
          spectrum_counts_(kEdepSpectrumBins, 0.0),
          spectrum_edep_(kEdepSpectrumBins, 0.0)
    {
      G4AccumulableManager* accumulable_manager = G4AccumulableManager::Instance();
      accumulable_manager->SetVerboseLevel(0);
      accumulable_manager->Register(edep_);
      accumulable_manager->Register(spectrum_underflow_);
      accumulable_manager->Register(spectrum_overflow_);
      for (auto& component_edep : component_edep_) {
        component_edep = std::make_unique<G4Accumulable<G4double>>(0.0);
        accumulable_manager->Register(*component_edep);
      }
      accumulable_manager->Register(spectrum_counts_);
      accumulable_manager->Register(spectrum_edep_);
    }

    void BeginOfRunAction(const G4Run*) override
    {
      // reset run accumulators
      G4RunManager::GetRunManager()->SetRandomNumberStore(false);
      G4AccumulableManager::Instance()->Reset();
      edep_spectrum_edges_mev_.resize(kEdepSpectrumBins + 1);
      edep_spectrum_edges_mev_[0] = 0.0;
      for (std::size_t i = 1; i <= kEdepSpectrumBins; ++i) {
        const G4double fraction =
          static_cast<G4double>(i - 1) / static_cast<G4double>(kEdepSpectrumBins - 1);
        edep_spectrum_edges_mev_[i] =
          kEdepSpectrumMinMeV *
          std::pow(kEdepSpectrumMaxMeV / kEdepSpectrumMinMeV, fraction);
      }
      if (!record_results_) {
        return;
      }
      results_.edep_spectrum_edges_mev.assign(
        edep_spectrum_edges_mev_.begin(), edep_spectrum_edges_mev_.end());
      results_.edep_spectrum_counts.assign(kEdepSpectrumBins, 0);
      results_.edep_spectrum_edep_mev.assign(kEdepSpectrumBins, 0.0);
      results_.edep_spectrum_underflow = 0;
      results_.edep_spectrum_overflow = 0;
      results_.component_names.clear();
      results_.component_edep_mev.assign(component_edep_.size(), 0.0);
      results_.component_mass_kg.assign(component_edep_.size(), 0.0);
      results_.component_dose_gy.assign(component_edep_.size(), 0.0);
    }

    void EndOfRunAction(const G4Run* run) override
    {
      if (!record_results_) {
        return;
      }
      const G4int nof_events = run->GetNumberOfEvent();
      results_.last_events_run = static_cast<std::size_t>(nof_events);
      if (nof_events == 0) {
        results_.last_total_edep_mev = 0.0;
        results_.last_dose_gy = 0.0;
        return;
      }

      // collect dose from accumulated energy deposition
      G4AccumulableManager::Instance()->Merge();
      const G4double edep = edep_.GetValue();
      const auto* det_construction = static_cast<const DetectorConstruction*>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());
      const auto& scoring_volumes = det_construction->GetScoringVolumes();
      results_.component_names = det_construction->GetScoringNames();
      results_.component_edep_mev.assign(component_edep_.size(), 0.0);
      results_.component_mass_kg.assign(component_edep_.size(), 0.0);
      results_.component_dose_gy.assign(component_edep_.size(), 0.0);
      results_.edep_spectrum_counts.assign(kEdepSpectrumBins, 0);
      results_.edep_spectrum_edep_mev.assign(kEdepSpectrumBins, 0.0);
      results_.edep_spectrum_underflow =
        static_cast<std::size_t>(std::llround(spectrum_underflow_.GetValue()));
      results_.edep_spectrum_overflow =
        static_cast<std::size_t>(std::llround(spectrum_overflow_.GetValue()));
      for (std::size_t i = 0; i < kEdepSpectrumBins; ++i) {
        results_.edep_spectrum_counts[i] =
          static_cast<std::size_t>(std::llround(spectrum_counts_[i]));
        results_.edep_spectrum_edep_mev[i] = spectrum_edep_[i];
      }

      G4double scored_mass = 0.0;
      for (std::size_t i = 0; i < component_edep_.size(); ++i) {
        const G4double component_edep_value = component_edep_[i]->GetValue();
        const G4double component_mass =
          (i < scoring_volumes.size()) ? scoring_volumes[i]->GetMass(false, false) : 0.0;
        scored_mass += component_mass;
        results_.component_edep_mev[i] = component_edep_value / MeV;
        results_.component_mass_kg[i] = component_mass / kg;
        results_.component_dose_gy[i] =
          (component_mass > 0.0) ? (component_edep_value / component_mass) / gray : 0.0;
      }
      const G4double dose = (scored_mass > 0.0) ? (edep / scored_mass) : 0.0;

      results_.last_total_edep_mev = edep / MeV;
      results_.last_dose_gy = dose / gray;
    }

    void AddEventEdep(
      G4double raw_edep,
      G4double weighted_edep,
      const std::vector<G4double>& component_weighted_edep)
    {
      edep_ += weighted_edep;
      for (std::size_t i = 0; i < component_weighted_edep.size(); ++i) {
        *component_edep_[i] += component_weighted_edep[i];
      }

      // bin physical event edep, and store weighted edep as the scored value
      const G4double raw_edep_mev = raw_edep / MeV;
      if (raw_edep_mev < 0.0) {
        spectrum_underflow_ += 1.0;
        return;
      }
      if (raw_edep_mev >= kEdepSpectrumMaxMeV) {
        spectrum_overflow_ += 1.0;
        return;
      }

      const auto upper_edge = std::upper_bound(
        edep_spectrum_edges_mev_.begin(), edep_spectrum_edges_mev_.end(), raw_edep_mev);
      const auto bin = static_cast<std::size_t>(
        std::distance(edep_spectrum_edges_mev_.begin(), upper_edge) - 1);
      spectrum_counts_[bin] += 1.0;
      spectrum_edep_[bin] += weighted_edep / MeV;
    }

  private:
    Results& results_;
    bool record_results_;
    G4Accumulable<G4double> edep_ = 0.0;
    G4Accumulable<G4double> spectrum_underflow_ = 0.0;
    G4Accumulable<G4double> spectrum_overflow_ = 0.0;
    std::vector<std::unique_ptr<G4Accumulable<G4double>>> component_edep_;
    G4AccVector<G4double> spectrum_counts_;
    G4AccVector<G4double> spectrum_edep_;
    std::vector<G4double> edep_spectrum_edges_mev_;
};

class BridgeEventAction : public G4UserEventAction
{
  public:
    BridgeEventAction(BridgeRunAction& run_action, std::size_t n_components)
        : run_action_(run_action),
          component_raw_edep_(n_components, 0.0),
          component_weighted_edep_(n_components, 0.0)
    {
    }

    void BeginOfEventAction(const G4Event*) override
    {
      // reset event energy deposition
      raw_edep_ = 0.0;
      edep_ = 0.0;
      std::fill(component_raw_edep_.begin(), component_raw_edep_.end(), 0.0);
      std::fill(component_weighted_edep_.begin(), component_weighted_edep_.end(), 0.0);
    }

    void EndOfEventAction(const G4Event*) override
    {
      // send event energy deposition to the run accumulator
      run_action_.AddEventEdep(raw_edep_, edep_, component_weighted_edep_);
    }

    void AddEdep(std::size_t component_idx, G4double raw_edep, G4double weighted_edep)
    {
      raw_edep_ += raw_edep;
      edep_ += weighted_edep;
      component_raw_edep_[component_idx] += raw_edep;
      component_weighted_edep_[component_idx] += weighted_edep;
    }

  private:
    BridgeRunAction& run_action_;
    G4double raw_edep_ = 0.0;
    G4double edep_ = 0.0;
    std::vector<G4double> component_raw_edep_;
    std::vector<G4double> component_weighted_edep_;
};

class BridgeSteppingAction : public G4UserSteppingAction
{
  public:
    explicit BridgeSteppingAction(BridgeEventAction& event_action) : event_action_(event_action) {}

    void UserSteppingAction(const G4Step* step) override
    {
      // cache scoring volumes on first step
      if (!scoring_volumes_cached_) {
        const auto* det_construction = static_cast<const DetectorConstruction*>(
          G4RunManager::GetRunManager()->GetUserDetectorConstruction());
        const auto& scoring_volumes = det_construction->GetScoringVolumes();
        for (std::size_t i = 0; i < scoring_volumes.size(); ++i) {
          scoring_volume_indices_[scoring_volumes[i]] = i;
        }
        scoring_volumes_cached_ = true;
      }

      G4LogicalVolume* volume =
        step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
      const auto component = scoring_volume_indices_.find(volume);
      if (component == scoring_volume_indices_.end()) {
        return;
      }

      // score weighted energy deposition inside the component
      const G4double edep_step = step->GetTotalEnergyDeposit();
      const G4double weighted_edep_step = edep_step * step->GetPreStepPoint()->GetWeight();
      event_action_.AddEdep(component->second, edep_step, weighted_edep_step);
    }

  private:
    BridgeEventAction& event_action_;
    std::unordered_map<G4LogicalVolume*, std::size_t> scoring_volume_indices_;
    bool scoring_volumes_cached_ = false;
};

class BridgeActionInitialization : public G4VUserActionInitialization
{
  public:
    BridgeActionInitialization(
      const PrimaryBank& primary_bank,
      const SourceDistribution& source_distribution,
      Results& results,
      std::size_t n_scoring_components,
      bool multithreaded)
        : primary_bank_(primary_bank),
          source_distribution_(source_distribution),
          results_(results),
          n_scoring_components_(n_scoring_components),
          multithreaded_(multithreaded)
    {
    }

    void BuildForMaster() const override
    {
      SetUserAction(new BridgeRunAction(results_, n_scoring_components_, true));
    }

    void Build() const override
    {
      // install bridge actions for worker runs
      SetUserAction(new BankPrimaryGeneratorAction(primary_bank_, source_distribution_));
      auto* run_action =
        new BridgeRunAction(results_, n_scoring_components_, !multithreaded_);
      SetUserAction(run_action);
      auto* event_action = new BridgeEventAction(*run_action, n_scoring_components_);
      SetUserAction(event_action);
      SetUserAction(new BridgeSteppingAction(*event_action));
    }

  private:
    const PrimaryBank& primary_bank_;
    const SourceDistribution& source_distribution_;
    Results& results_;
    std::size_t n_scoring_components_;
    bool multithreaded_;
};

class SilentUIsession : public G4UIsession
{
  public:
    G4int ReceiveG4cout(const G4String&) override { return 0; }
    G4int ReceiveG4cerr(const G4String&) override { return 0; }
};

}  // namespace

Session::Session(SessionConfig config)
    : config_(std::move(config)),
      physics_list_name_(ResolvePhysicsList(config_)),
      primary_bank_(std::make_unique<PrimaryBank>()),
      source_distribution_(std::make_unique<SourceDistribution>())
{
  ValidateGeometryConfig(config_);
  results_.physics_list = physics_list_name_;
}

Session::~Session()
{
  try {
    close();
  }
  catch (...) {
    // Destructors must not throw.
  }
}

void Session::initialize()
{
  if (initialized_) {
    return;
  }

  // create Geant4 run manager
  if (config_.n_threads > 1) {
#ifdef G4MULTITHREADED
    run_manager_.reset(G4RunManagerFactory::CreateRunManager(G4RunManagerType::MTOnly));
    auto* mt_run_manager = dynamic_cast<G4MTRunManager*>(run_manager_.get());
    if (!mt_run_manager) {
      throw std::runtime_error("Failed to create Geant4 MT run manager.");
    }
    mt_run_manager->SetNumberOfThreads(config_.n_threads);
#else
    throw std::runtime_error(
      "Geant4 bridge was built without multithreading support; n_threads > 1 requested.");
#endif
  } else {
    run_manager_.reset(G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly));
  }
  if (!run_manager_) {
    throw std::runtime_error("Failed to create Geant4 run manager.");
  }

  // seed the worker before Geant4 builds or samples anything
  CLHEP::HepRandom::setTheSeed(config_.random_seed);

  // configure detector geometry
  run_manager_->SetUserInitialization(
    new DetectorConstruction(
      config_.world_size_mm,
      config_.detector_size_mm,
      config_.envelope_material,
      config_.device_components));

  // configure physics list
  G4PhysListFactory phys_factory;
  G4VModularPhysicsList* physics_list =
    phys_factory.GetReferencePhysList(physics_list_name_);
  if (!physics_list) {
    std::ostringstream msg;
    msg << "Invalid physics_list '" << physics_list_name_
        << "'. Could not create reference Geant4 physics list.";
    throw std::runtime_error(msg.str());
  }
  physics_list->SetVerboseLevel(0);
  G4EmParameters::Instance()->SetVerbose(0);
  G4HadronicParameters::Instance()->SetVerboseLevel(0);
  run_manager_->SetUserInitialization(physics_list);

  // install bridge action stack
  run_manager_->SetUserInitialization(new BridgeActionInitialization(
    *primary_bank_,
    *source_distribution_,
    results_,
    CountScoringComponents(config_.device_components),
    config_.n_threads > 1));
  run_manager_->Initialize();

  // silence Geant4 command output
  G4UImanager::GetUIpointer()->ApplyCommand("/control/verbose 0");
  G4UImanager::GetUIpointer()->ApplyCommand("/run/verbose 0");
  G4UImanager::GetUIpointer()->ApplyCommand("/event/verbose 0");
  G4UImanager::GetUIpointer()->ApplyCommand("/tracking/verbose 0");

  initialized_ = true;
  results_.initialized = true;
  results_.status = "initialized";
}

void Session::load_primaries(
  const py::array_t<double, py::array::c_style | py::array::forcecast>& primaries)
{
  if (!initialized_) {
    throw std::runtime_error("Session must be initialized before load_primaries().");
  }

  primary_bank_->Load(primaries);
  source_distribution_->Clear();

  // update active source summary
  has_source_ = true;
  results_.has_source = true;
  results_.loaded_primaries = primary_bank_->Size();
  results_.status = "primaries_loaded";
}

void Session::load_source_distribution(
  const py::array_t<double, py::array::c_style | py::array::forcecast>& box_bounds_mm,
  const py::array_t<double, py::array::c_style | py::array::forcecast>& mu_edges,
  const py::array_t<double, py::array::c_style | py::array::forcecast>& azi_edges,
  const py::array_t<double, py::array::c_style | py::array::forcecast>& energy_edges_mev,
  const py::array_t<double, py::array::c_style | py::array::forcecast>& weights,
  std::size_t n_u,
  std::size_t n_v,
  std::size_t n_events)
{
  if (!initialized_) {
    throw std::runtime_error(
      "Session must be initialized before load_source_distribution().");
  }

  // load sampled distribution source instead of exact bank rows
  source_distribution_->Load(
    box_bounds_mm, mu_edges, azi_edges, energy_edges_mev, weights, n_u, n_v, n_events);
  primary_bank_->Clear();

  // update active source summary
  has_source_ = true;
  results_.has_source = true;
  results_.loaded_primaries = source_distribution_->Size();
  results_.status = "distribution_loaded";
}

void Session::beam_on()
{
  if (!initialized_) {
    throw std::runtime_error("Session must be initialized before beam_on().");
  }
  if (!has_source_) {
    throw std::runtime_error(
      "No source loaded. Call load_primaries() or load_source_distribution() before beam_on().");
  }

  // choose event count from the active source
  const std::size_t n_events =
    source_distribution_->Loaded() ? source_distribution_->Size() : primary_bank_->Size();
  if (n_events > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::runtime_error("Geant4 source too large for BeamOn(int) event count.");
  }

  // run Geant4 events
  run_manager_->BeamOn(static_cast<int>(n_events));
  results_.status = "ok";

  std::cout << "[geant4_bridge] beam_on loaded_rows=" << results_.loaded_primaries
            << " events_run=" << results_.last_events_run
            << " physics_list=" << results_.physics_list
            << " status=" << results_.status << "\n";
}

Results Session::get_results() const
{
  return results_;
}

void Session::close()
{
  if (!run_manager_) {
    return;
  }

  // silence teardown output before destroying the run manager
  G4UImanager* ui = G4UImanager::GetUIpointer();
  if (ui) {
    if (!silent_ui_session_) {
      silent_ui_session_ = std::make_unique<SilentUIsession>();
    }
    ui->SetCoutDestination(silent_ui_session_.get());
  }

  // reset bridge state
  run_manager_.reset();
  initialized_ = false;
  has_source_ = false;
  primary_bank_->Clear();
  source_distribution_->Clear();
  results_.ResetRuntime();
  results_.initialized = false;
  results_.has_source = false;
  results_.loaded_primaries = 0;
  results_.physics_list = physics_list_name_;
  results_.status = "closed";
}

void Session::ValidateGeometryConfig(const SessionConfig& config)
{
  // validate geometry vector entries
  auto validate_vec = [](const std::array<double, 3>& v, const char* name) {
    for (int i = 0; i < 3; ++i) {
      if (!std::isfinite(v[static_cast<std::size_t>(i)])) {
        std::ostringstream msg;
        msg << "config." << name << " has non-finite value at index " << i << ".";
        throw std::runtime_error(msg.str());
      }
      if (v[static_cast<std::size_t>(i)] <= 0.0) {
        std::ostringstream msg;
        msg << "config." << name << " must have strictly positive entries.";
        throw std::runtime_error(msg.str());
      }
    }
  };
  validate_vec(config.world_size_mm, "world_size_mm");
  validate_vec(config.detector_size_mm, "detector_size_mm");
  if (config.detector_material.empty()) {
    throw std::runtime_error("config.detector_material cannot be empty.");
  }
  if (config.envelope_material.empty()) {
    throw std::runtime_error("config.envelope_material cannot be empty.");
  }
  if (config.device_components.empty()) {
    throw std::runtime_error("config.device_components cannot be empty.");
  }
  if (config.random_seed <= 0) {
    throw std::runtime_error("config.random_seed must be positive.");
  }
  if (config.n_threads < 1) {
    throw std::runtime_error("config.n_threads must be at least 1.");
  }

  // require detector to fit inside world
  for (int i = 0; i < 3; ++i) {
    const auto idx = static_cast<std::size_t>(i);
    if (config.detector_size_mm[idx] >= config.world_size_mm[idx]) {
      std::ostringstream msg;
      msg << "config.detector_size_mm[" << i
          << "] must be smaller than config.world_size_mm[" << i << "].";
      throw std::runtime_error(msg.str());
    }
  }

  std::unordered_map<std::string, std::size_t> component_indices;
  component_indices.reserve(config.device_components.size());
  for (std::size_t component_idx = 0; component_idx < config.device_components.size();
       ++component_idx) {
    const auto& component = config.device_components[component_idx];
    if (component.name.empty()) {
      throw std::runtime_error("Device component name cannot be empty.");
    }
    if (!component_indices.emplace(component.name, component_idx).second) {
      std::ostringstream msg;
      msg << "Duplicate device component name '" << component.name << "'.";
      throw std::runtime_error(msg.str());
    }
    if (component.material.empty()) {
      std::ostringstream msg;
      msg << "Device component '" << component.name << "' material cannot be empty.";
      throw std::runtime_error(msg.str());
    }
    for (int axis = 0; axis < 3; ++axis) {
      const auto idx = static_cast<std::size_t>(axis);
      const double center = component.center_mm[idx];
      const double size = component.size_mm[idx];
      if (!std::isfinite(center) || !std::isfinite(size)) {
        std::ostringstream msg;
        msg << "Device component '" << component.name << "' has non-finite geometry.";
        throw std::runtime_error(msg.str());
      }
      if (size <= 0.0) {
        std::ostringstream msg;
        msg << "Device component '" << component.name
            << "' size_mm entries must be positive.";
        throw std::runtime_error(msg.str());
      }
    }
  }

  for (const auto& component : config.device_components) {
    if (component.parent == component.name) {
      std::ostringstream msg;
      msg << "Device component '" << component.name << "' cannot parent itself.";
      throw std::runtime_error(msg.str());
    }
    if (!component.parent.empty() && component_indices.count(component.parent) == 0) {
      std::ostringstream msg;
      msg << "Device component '" << component.name << "' references unknown parent '"
          << component.parent << "'.";
      throw std::runtime_error(msg.str());
    }
  }

  for (const auto& component : config.device_components) {
    std::unordered_set<std::string> seen;
    std::string parent = component.parent;
    while (!parent.empty()) {
      if (!seen.insert(parent).second) {
        std::ostringstream msg;
        msg << "Device component parent cycle includes '" << parent << "'.";
        throw std::runtime_error(msg.str());
      }
      parent = config.device_components[component_indices.at(parent)].parent;
    }
  }

  for (const auto& component : config.device_components) {
    const DeviceComponent* parent = nullptr;
    if (!component.parent.empty()) {
      parent = &config.device_components[component_indices.at(component.parent)];
    }
    for (int axis = 0; axis < 3; ++axis) {
      const auto idx = static_cast<std::size_t>(axis);
      const double parent_center = parent ? parent->center_mm[idx] : 0.0;
      const double parent_size = parent ? parent->size_mm[idx] : config.detector_size_mm[idx];
      const double distance = std::abs(component.center_mm[idx] - parent_center);
      const bool outside_parent =
        distance + 0.5 * component.size_mm[idx] >
        0.5 * parent_size + kGeometryToleranceMm;
      if (outside_parent) {
        std::ostringstream msg;
        msg << "Device component '" << component.name << "' does not fit inside "
            << (parent ? "parent" : "detector envelope") << ".";
        throw std::runtime_error(msg.str());
      }
    }
  }

  for (std::size_t i = 0; i < config.device_components.size(); ++i) {
    for (std::size_t j = i + 1; j < config.device_components.size(); ++j) {
      const auto& a = config.device_components[i];
      const auto& b = config.device_components[j];
      if (a.parent != b.parent) {
        continue;
      }
      bool overlap = true;
      for (int axis = 0; axis < 3; ++axis) {
        const auto idx = static_cast<std::size_t>(axis);
        const double distance = std::abs(a.center_mm[idx] - b.center_mm[idx]);
        const double half_sum = 0.5 * (a.size_mm[idx] + b.size_mm[idx]);
        if (distance + kGeometryToleranceMm >= half_sum) {
          overlap = false;
          break;
        }
      }
      if (overlap) {
        std::ostringstream msg;
        msg << "Device components '" << a.name << "' and '" << b.name
            << "' overlap.";
        throw std::runtime_error(msg.str());
      }
    }
  }
}

std::string Session::ResolvePhysicsList(const SessionConfig& config) const
{
  // validate physics list name
  const std::string name = config.physics_list;
  if (name.empty()) {
    throw std::runtime_error("config.physics_list cannot be empty.");
  }
#if defined(__APPLE__)
  if (name.find("_HP") != std::string::npos) {
    std::cout
      << "[geant4_bridge] warning: macOS + Geant4 *_HP physics lists may emit non-critical "
         "G4Cache mutex-lock teardown warnings due to an upstream Geant4 issue.\n";
  }
#endif
  return name;
}

}  // namespace g4bridge
