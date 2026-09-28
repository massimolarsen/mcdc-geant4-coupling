#include "session.hpp"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include "bank_primary_generator_action.hpp"
#include "detector_construction.hpp"
#include "primary_bank.hpp"
#include "seu_scoring.hpp"
#include "source_distribution.hpp"

#include "G4Accumulable.hh"
#include "G4AccumulableManager.hh"
#include "G4AccVector.hh"
#include "G4EmParameters.hh"
#include "G4Event.hh"
#include "G4HadronicParameters.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4MaterialCutsCouple.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#endif
#include "G4PhysListFactory.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4ProductionCuts.hh"
#include "G4ProductionCutsTable.hh"
#include "G4Region.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4Step.hh"
#include "G4Threading.hh"
#include "G4Track.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "G4UIsession.hh"
#include "G4UserEventAction.hh"
#include "G4UserRunAction.hh"
#include "G4UserSteppingAction.hh"
#include "G4UserTrackingAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VProcess.hh"
#include "G4Version.hh"
#include "CLHEP/Random/Random.h"

namespace py = pybind11;

namespace g4bridge
{

struct ProgressState
{
  std::atomic<std::size_t> completed{0};
  std::size_t total = 0;
  std::chrono::steady_clock::time_point started;
  std::mutex output_mutex;
};

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
          spectrum_edep_(kEdepSpectrumBins, 0.0),
          component_niel_(n_components, 0.0),
          component_ionizing_(n_components, 0.0),
          species_ionizing_(n_components * kSeUSpeciesCount, 0.0),
          primary_ionizing_(n_components, 0.0),
          secondary_ionizing_(n_components, 0.0),
          ionizing_count_(n_components * kSeUStride, 0.0),
          ionizing_sumw_(n_components * kSeUStride, 0.0),
          ionizing_sumw2_(n_components * kSeUStride, 0.0)
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
      accumulable_manager->Register(component_niel_);
      accumulable_manager->Register(component_ionizing_);
      accumulable_manager->Register(species_ionizing_);
      accumulable_manager->Register(primary_ionizing_);
      accumulable_manager->Register(secondary_ionizing_);
      accumulable_manager->Register(ionizing_count_);
      accumulable_manager->Register(ionizing_sumw_);
      accumulable_manager->Register(ionizing_sumw2_);
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
      ionizing_edges_mev_ = SeUEdgesMeV();
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
      results_.component_niel_mev.assign(component_edep_.size(), 0.0);
      results_.component_ionizing_mev.assign(component_edep_.size(), 0.0);
      results_.seu_species_names = SeUSpeciesNames();
      results_.component_species_ionizing_mev.assign(species_ionizing_.size(), 0.0);
      results_.component_primary_ionizing_mev.assign(component_edep_.size(), 0.0);
      results_.component_secondary_ionizing_mev.assign(component_edep_.size(), 0.0);
      results_.component_event_ionizing_edges_mev = ionizing_edges_mev_;
      results_.component_event_ionizing_count.assign(ionizing_count_.size(), 0.0);
      results_.component_event_ionizing_sumw.assign(ionizing_sumw_.size(), 0.0);
      results_.component_event_ionizing_sumw2.assign(ionizing_sumw2_.size(), 0.0);
    }

    void EndOfRunAction(const G4Run* run) override
    {
      G4AccumulableManager::Instance()->Merge();
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
      const G4double edep = edep_.GetValue();
      const auto* det_construction = static_cast<const DetectorConstruction*>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());
      const auto& scoring_volumes = det_construction->GetScoringVolumes();
      results_.component_names = det_construction->GetScoringNames();
      results_.component_edep_mev.assign(component_edep_.size(), 0.0);
      results_.component_mass_kg.assign(component_edep_.size(), 0.0);
      results_.component_dose_gy.assign(component_edep_.size(), 0.0);
      results_.component_niel_mev.assign(component_edep_.size(), 0.0);
      results_.component_ionizing_mev.assign(component_edep_.size(), 0.0);
      results_.component_primary_ionizing_mev.assign(component_edep_.size(), 0.0);
      results_.component_secondary_ionizing_mev.assign(component_edep_.size(), 0.0);
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
        const G4double component_mass = scoring_volumes[i]->GetMass(false, false);
        scored_mass += component_mass;
        results_.component_edep_mev[i] = component_edep_value / MeV;
        results_.component_mass_kg[i] = component_mass / kg;
        results_.component_dose_gy[i] =
          (component_mass > 0.0) ? (component_edep_value / component_mass) / gray : 0.0;
        results_.component_niel_mev[i] = component_niel_[i] / MeV;
        results_.component_ionizing_mev[i] = component_ionizing_[i] / MeV;
        results_.component_primary_ionizing_mev[i] = primary_ionizing_[i] / MeV;
        results_.component_secondary_ionizing_mev[i] = secondary_ionizing_[i] / MeV;
      }
      results_.seu_species_names = SeUSpeciesNames();
      for (std::size_t i = 0; i < species_ionizing_.size(); ++i) {
        results_.component_species_ionizing_mev[i] = species_ionizing_[i] / MeV;
      }
      for (std::size_t i = 0; i < ionizing_count_.size(); ++i) {
        results_.component_event_ionizing_count[i] = ionizing_count_[i];
        results_.component_event_ionizing_sumw[i] = ionizing_sumw_[i];
        results_.component_event_ionizing_sumw2[i] = ionizing_sumw2_[i];
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

    void AddSeUEvent(
      const std::vector<G4double>& total,
      const std::vector<G4double>& niel,
      const std::vector<G4double>& species_ionizing,
      const std::vector<G4double>& primary_ionizing,
      const std::vector<G4double>& secondary_ionizing,
      G4double event_weight)
    {
      for (std::size_t i = 0; i < total.size(); ++i) {
        const G4double ion = total[i] - niel[i];
        component_niel_[i] += event_weight * niel[i];
        component_ionizing_[i] += event_weight * ion;
        primary_ionizing_[i] += event_weight * primary_ionizing[i];
        secondary_ionizing_[i] += event_weight * secondary_ionizing[i];
        const auto slot = i * kSeUStride + SeUBin(ion / MeV, ionizing_edges_mev_);
        ionizing_count_[slot] += 1.0;
        ionizing_sumw_[slot] += event_weight;
        ionizing_sumw2_[slot] += event_weight * event_weight;
      }
      for (std::size_t i = 0; i < species_ionizing.size(); ++i) {
        species_ionizing_[i] += event_weight * species_ionizing[i];
      }
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
    G4AccVector<G4double> component_niel_;
    G4AccVector<G4double> component_ionizing_;
    G4AccVector<G4double> species_ionizing_;
    G4AccVector<G4double> primary_ionizing_;
    G4AccVector<G4double> secondary_ionizing_;
    G4AccVector<G4double> ionizing_count_;
    G4AccVector<G4double> ionizing_sumw_;
    G4AccVector<G4double> ionizing_sumw2_;
    std::vector<G4double> edep_spectrum_edges_mev_;
    std::vector<double> ionizing_edges_mev_;
};

class BridgeEventAction : public G4UserEventAction
{
  public:
    BridgeEventAction(
      BridgeRunAction& run_action,
      std::size_t n_components,
      bool record_details,
      double min_ionizing_mev,
      std::string diagnostic_dir,
      std::shared_ptr<ProgressState> progress_state)
        : run_action_(run_action),
          component_weighted_edep_(n_components, 0.0),
          total_(n_components, 0.0),
          niel_(n_components, 0.0),
          species_ionizing_(n_components * kSeUSpeciesCount, 0.0),
          primary_ionizing_(n_components, 0.0),
          secondary_ionizing_(n_components, 0.0),
          record_details_(record_details),
          min_ionizing_mev_(min_ionizing_mev),
          diagnostic_dir_(std::move(diagnostic_dir)),
          progress_state_(std::move(progress_state))
    {
    }

    void BeginOfEventAction(const G4Event*) override
    {
      // reset event energy deposition
      raw_edep_ = 0.0;
      edep_ = 0.0;
      std::fill(component_weighted_edep_.begin(), component_weighted_edep_.end(), 0.0);
      std::fill(total_.begin(), total_.end(), 0.0);
      std::fill(niel_.begin(), niel_.end(), 0.0);
      std::fill(species_ionizing_.begin(), species_ionizing_.end(), 0.0);
      std::fill(primary_ionizing_.begin(), primary_ionizing_.end(), 0.0);
      std::fill(secondary_ionizing_.begin(), secondary_ionizing_.end(), 0.0);
      em_count_.fill(0);
      em_energy_mev_.fill(0.0);
      entry_rows_.clear();
      nuclear_rows_.clear();
    }

    void EndOfEventAction(const G4Event* event) override
    {
      // send event energy deposition to the run accumulator
      run_action_.AddEventEdep(raw_edep_, edep_, component_weighted_edep_);
      const G4double weight = event->GetPrimaryVertex(0)->GetPrimary()->GetWeight();
      run_action_.AddSeUEvent(
        total_, niel_, species_ionizing_, primary_ionizing_, secondary_ionizing_, weight);

      const auto completed = progress_state_->completed.fetch_add(1) + 1;
      const auto percent = 100 * completed / progress_state_->total;
      if (percent > 100 * (completed - 1) / progress_state_->total) {
        const auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(
          std::chrono::steady_clock::now() - progress_state_->started).count();
        std::lock_guard<std::mutex> lock(progress_state_->output_mutex);
        std::cerr << "Geant4 progress: " << completed << "/" << progress_state_->total
                  << " events (" << percent << "%) elapsed=" << elapsed << "s" << std::endl;
      }

      if (!record_details_) return;
      G4double largest = 0.0;
      for (std::size_t i = 0; i < total_.size(); ++i) {
        largest = std::max(largest, total_[i] - niel_[i]);
      }
      if (largest / MeV < min_ionizing_mev_) return;

      if (!diagnostics_.is_open()) {
        const auto thread_id = G4Threading::G4GetThreadId();
        diagnostics_.open(
          diagnostic_dir_ + "/worker_" + std::to_string(thread_id) + ".tsv",
          std::ios::app);
        if (!diagnostics_) throw std::runtime_error("Cannot open Geant4 SEU diagnostic stream.");
        diagnostics_ << std::setprecision(17);
      }
      const auto* run = G4RunManager::GetRunManager()->GetCurrentRun();
      const int run_id = run->GetRunID();
      const int event_id = event->GetEventID();
      for (std::size_t i = 0; i < total_.size(); ++i) {
        diagnostics_ << "E\t" << run_id << '\t' << event_id << '\t' << i << '\t'
                     << weight << '\t' << total_[i] / MeV << '\t' << niel_[i] / MeV
                     << '\t' << (total_[i] - niel_[i]) / MeV << '\t'
                     << primary_ionizing_[i] / MeV << '\t' << secondary_ionizing_[i] / MeV;
        for (std::size_t s = 0; s < kSeUSpeciesCount; ++s) {
          diagnostics_ << '\t' << species_ionizing_[i * kSeUSpeciesCount + s] / MeV;
        }
        diagnostics_ << '\n';
      }
      diagnostics_ << "M\t" << run_id << '\t' << event_id;
      for (std::size_t i = 0; i < 3; ++i) {
        diagnostics_ << '\t' << em_count_[i] << '\t' << em_energy_mev_[i];
      }
      diagnostics_ << '\n';
      for (const auto& row : entry_rows_) {
        diagnostics_ << "I\t" << run_id << '\t' << event_id << '\t' << row << '\n';
      }
      for (const auto& row : nuclear_rows_) {
        diagnostics_ << "N\t" << run_id << '\t' << event_id << '\t' << row << '\n';
      }
    }

    void AddEdep(std::size_t component_idx, G4double raw_edep, G4double weighted_edep)
    {
      raw_edep_ += raw_edep;
      edep_ += weighted_edep;
      component_weighted_edep_[component_idx] += weighted_edep;
    }

    void AddSeUStep(std::size_t component_idx, const G4Step* step)
    {
      const G4double total = step->GetTotalEnergyDeposit();
      const G4double niel = step->GetNonIonizingEnergyDeposit();
      if (!std::isfinite(total) || !std::isfinite(niel) || total < -1.e-12 * MeV ||
          niel < -1.e-12 * MeV || niel > total + 1.e-12 * MeV) {
        throw std::runtime_error("Geant4 step has invalid total or non-ionizing deposition.");
      }
      const G4double ion = total - niel;
      total_[component_idx] += total;
      niel_[component_idx] += niel;
      const G4Track* track = step->GetTrack();
      const int pdg = track->GetDefinition()->GetPDGEncoding();
      species_ionizing_[component_idx * kSeUSpeciesCount + SeUSpeciesForPdg(pdg)] += ion;
      if (track->GetParentID() == 0) {
        primary_ionizing_[component_idx] += ion;
      } else {
        secondary_ionizing_[component_idx] += ion;
      }
    }

    void AddEntry(std::size_t component_idx, const G4Step* step)
    {
      if (!record_details_) return;
      const G4Track* track = step->GetTrack();
      const auto* point = step->GetPreStepPoint();
      const auto& position = point->GetPosition();
      const auto& direction = point->GetMomentumDirection();
      const auto& vertex = track->GetVertexPosition();
      const auto* origin = track->GetLogicalVolumeAtVertex();
      const auto* creator = track->GetCreatorProcess();
      std::ostringstream row;
      row << std::setprecision(17) << component_idx << '\t' << track->GetTrackID() << '\t'
          << track->GetParentID() << '\t' << track->GetDefinition()->GetPDGEncoding() << '\t'
          << point->GetKineticEnergy() / MeV << '\t' << position.x() / mm << '\t'
          << position.y() / mm << '\t' << position.z() / mm << '\t'
          << direction.x() << '\t' << direction.y() << '\t' << direction.z() << '\t'
          << (creator ? creator->GetProcessName() : "primary") << '\t'
          << vertex.x() / mm << '\t' << vertex.y() / mm << '\t' << vertex.z() / mm
          << '\t' << track->GetVertexKineticEnergy() / MeV << '\t'
          << (origin ? origin->GetName() : "");
      entry_rows_.push_back(row.str());
    }

    void AddSecondary(const G4Track* track)
    {
      if (!record_details_ || track->GetParentID() == 0) return;
      const int pdg = track->GetDefinition()->GetPDGEncoding();
      const int code = std::abs(pdg);
      const double energy_mev = track->GetVertexKineticEnergy() / MeV;
      if (code == 11 || code == 22) {
        const std::size_t i = code == 22 ? 2 : (pdg == 11 ? 0 : 1);
        ++em_count_[i];
        em_energy_mev_[i] += energy_mev;
        return;
      }
      if (code != 2112 && code != 2212 && code < 1000000000) return;
      const auto& vertex = track->GetVertexPosition();
      const auto* origin = track->GetLogicalVolumeAtVertex();
      const auto* creator = track->GetCreatorProcess();
      std::ostringstream row;
      row << std::setprecision(17) << track->GetTrackID() << '\t' << track->GetParentID()
          << '\t' << pdg << '\t' << energy_mev << '\t' << vertex.x() / mm << '\t'
          << vertex.y() / mm << '\t' << vertex.z() / mm << '\t'
          << (origin ? origin->GetName() : "") << '\t'
          << (creator ? creator->GetProcessName() : "");
      nuclear_rows_.push_back(row.str());
    }

  private:
    BridgeRunAction& run_action_;
    G4double raw_edep_ = 0.0;
    G4double edep_ = 0.0;
    std::vector<G4double> component_weighted_edep_;
    std::vector<G4double> total_;
    std::vector<G4double> niel_;
    std::vector<G4double> species_ionizing_;
    std::vector<G4double> primary_ionizing_;
    std::vector<G4double> secondary_ionizing_;
    bool record_details_;
    double min_ionizing_mev_;
    std::string diagnostic_dir_;
    std::shared_ptr<ProgressState> progress_state_;
    std::ofstream diagnostics_;
    std::array<std::size_t, 3> em_count_{};
    std::array<double, 3> em_energy_mev_{};
    std::vector<std::string> entry_rows_;
    std::vector<std::string> nuclear_rows_;
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
        const auto& scoring_names = det_construction->GetScoringNames();
        for (std::size_t i = 0; i < scoring_names.size(); ++i) {
          scoring_volume_indices_[scoring_names[i]] = i;
        }
        scoring_volumes_cached_ = true;
      }

      G4LogicalVolume* volume =
        step->GetPreStepPoint()->GetTouchableHandle()->GetVolume()->GetLogicalVolume();
      // In Geant4 MT runs, worker geometry can use cloned logical-volume
      // pointers, so match by stable component name instead of pointer value.
      const auto component = scoring_volume_indices_.find(volume->GetName());
      if (component == scoring_volume_indices_.end()) {
        return;
      }

      if (step->GetPreStepPoint()->GetStepStatus() == fGeomBoundary ||
          (step->GetTrack()->GetCurrentStepNumber() == 1 &&
           step->GetTrack()->GetParentID() == 0)) {
        event_action_.AddEntry(component->second, step);
      }
      event_action_.AddSeUStep(component->second, step);
      // score weighted energy deposition inside the component
      const G4double edep_step = step->GetTotalEnergyDeposit();
      const G4double weighted_edep_step = edep_step * step->GetPreStepPoint()->GetWeight();
      event_action_.AddEdep(component->second, edep_step, weighted_edep_step);
    }

  private:
    BridgeEventAction& event_action_;
    std::unordered_map<std::string, std::size_t> scoring_volume_indices_;
    bool scoring_volumes_cached_ = false;
};

class BridgeTrackingAction : public G4UserTrackingAction
{
  public:
    explicit BridgeTrackingAction(BridgeEventAction& event_action) : event_action_(event_action) {}

    void PreUserTrackingAction(const G4Track* track) override
    {
      event_action_.AddSecondary(track);
    }

  private:
    BridgeEventAction& event_action_;
};

class BridgeActionInitialization : public G4VUserActionInitialization
{
  public:
    BridgeActionInitialization(
      const PrimaryBank& primary_bank,
      const SourceDistribution& source_distribution,
      Results& results,
      std::size_t n_scoring_components,
      bool multithreaded,
      bool record_details,
      double min_ionizing_mev,
      std::string diagnostic_dir,
      std::shared_ptr<ProgressState> progress_state)
        : primary_bank_(primary_bank),
          source_distribution_(source_distribution),
          results_(results),
          n_scoring_components_(n_scoring_components),
          multithreaded_(multithreaded),
          record_details_(record_details),
          min_ionizing_mev_(min_ionizing_mev),
          diagnostic_dir_(std::move(diagnostic_dir)),
          progress_state_(std::move(progress_state))
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
      auto* event_action = new BridgeEventAction(
        *run_action, n_scoring_components_, record_details_, min_ionizing_mev_, diagnostic_dir_,
        progress_state_);
      SetUserAction(event_action);
      SetUserAction(new BridgeSteppingAction(*event_action));
      SetUserAction(new BridgeTrackingAction(*event_action));
    }

  private:
    const PrimaryBank& primary_bank_;
    const SourceDistribution& source_distribution_;
    Results& results_;
    std::size_t n_scoring_components_;
    bool multithreaded_;
    bool record_details_;
    double min_ionizing_mev_;
    std::string diagnostic_dir_;
    std::shared_ptr<ProgressState> progress_state_;
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
      source_distribution_(std::make_unique<SourceDistribution>()),
      progress_state_(std::make_shared<ProgressState>())
{
  ValidateGeometryConfig(config_);
  results_.physics_list = physics_list_name_;
  results_.geant4_version = G4Version;
  results_.em_production_cut_mm = config_.em_production_cut_mm;
  results_.proton_production_cut_mm = 0.0;
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
      config_.device_components,
      config_.em_production_cut_mm));

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
    config_.n_threads > 1,
    config_.record_seu_events,
    config_.diagnostic_min_Eion_mev,
    config_.diagnostic_dir,
    progress_state_));
  run_manager_->Initialize();

  // silence Geant4 command output
  G4UImanager::GetUIpointer()->ApplyCommand("/control/verbose 0");
  G4UImanager::GetUIpointer()->ApplyCommand("/run/verbose 0");
  G4UImanager::GetUIpointer()->ApplyCommand("/event/verbose 0");
  G4UImanager::GetUIpointer()->ApplyCommand("/tracking/verbose 0");

  initialized_ = true;
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
  results_.loaded_primaries = primary_bank_->Size();
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
  results_.loaded_primaries = source_distribution_->Size();
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
  progress_state_->completed = 0;
  progress_state_->total = n_events;
  progress_state_->started = std::chrono::steady_clock::now();
  run_manager_->BeamOn(static_cast<int>(n_events));

  // Geant4 populates the material-cuts table during the first BeamOn.
  results_.electronics_cut_materials.clear();
  results_.electronics_cut_energy_mev.clear();
  const auto* detector = static_cast<const DetectorConstruction*>(
    run_manager_->GetUserDetectorConstruction());
  if (const auto* electronics = detector->GetElectronicsRegion()) {
    const auto* table = G4ProductionCutsTable::GetProductionCutsTable();
    const std::array<const char*, 4> particles = {"gamma", "e-", "e+", "proton"};
    for (std::size_t i = 0; i < table->GetTableSize(); ++i) {
      const auto* couple = table->GetMaterialCutsCouple(static_cast<G4int>(i));
      if (couple->GetProductionCuts() != electronics->GetProductionCuts()) continue;
      results_.electronics_cut_materials.push_back(couple->GetMaterial()->GetName());
      for (const char* particle : particles) {
        const auto index = G4ProductionCuts::GetIndex(particle);
        results_.electronics_cut_energy_mev.push_back(
          (*table->GetEnergyCutsVector(static_cast<std::size_t>(index)))[i] / MeV);
      }
    }
  }
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
  results_.loaded_primaries = 0;
  results_.physics_list = physics_list_name_;
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
  if (!std::isfinite(config.em_production_cut_mm) || config.em_production_cut_mm < 0.0) {
    throw std::runtime_error("config.em_production_cut_mm must be finite and non-negative.");
  }
  if (!std::isfinite(config.diagnostic_min_Eion_mev) ||
      config.diagnostic_min_Eion_mev < 0.0) {
    throw std::runtime_error("config.diagnostic_min_Eion_mev must be finite and non-negative.");
  }
  if (config.record_seu_events && config.diagnostic_dir.empty()) {
    throw std::runtime_error("config.diagnostic_dir is required when recording SEU events.");
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
