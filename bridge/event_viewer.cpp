#include "event_viewer.hpp"

#include <map>
#include <memory>
#include <stdexcept>

#include "detector_construction.hpp"

#include "G4AttDef.hh"
#include "G4AttDefStore.hh"
#include "G4AttValue.hh"
#include "G4Electron.hh"
#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4Geantino.hh"
#include "G4ParticleGun.hh"
#include "G4Positron.hh"
#include "G4Proton.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4SystemOfUnits.hh"
#include "G4TrajectoryContainer.hh"
#include "G4UIcommand.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4UserEventAction.hh"
#include "G4VTrajectory.hh"
#include "G4VTrajectoryPoint.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VisExecutive.hh"

namespace g4bridge
{

namespace
{

class RecordedPoint : public G4VTrajectoryPoint
{
  public:
    RecordedPoint(const RecordedTracks& tracks, std::size_t i)
        : position_(tracks.x_mm[i] * mm, tracks.y_mm[i] * mm, tracks.z_mm[i] * mm),
          time_(tracks.t_ns[i] * ns),
          kinetic_energy_(tracks.ke_mev[i] * MeV),
          edep_(tracks.edep_mev[i] * MeV),
          process_(tracks.process[i]),
          volume_(tracks.volume[i])
    {
    }

    const G4ThreeVector GetPosition() const override { return position_; }

    const std::map<G4String, G4AttDef>* GetAttDefs() const override
    {
      G4bool is_new;
      auto* store = G4AttDefStore::GetInstance("RecordedPoint", is_new);
      if (is_new) {
        (*store)["Pos"] = G4AttDef("Pos", "Position", "Physics", "G4BestUnit", "G4ThreeVector");
        (*store)["T"] = G4AttDef("T", "Global time", "Physics", "G4BestUnit", "G4double");
        (*store)["KE"] = G4AttDef("KE", "Kinetic energy", "Physics", "G4BestUnit", "G4double");
        (*store)["EDep"] =
          G4AttDef("EDep", "Energy deposited in the step to here", "Physics", "G4BestUnit",
                   "G4double");
        (*store)["Proc"] = G4AttDef("Proc", "Process that ended the step", "Physics", "", "G4String");
        (*store)["Vol"] = G4AttDef("Vol", "Volume of the step", "Physics", "", "G4String");
      }
      return store;
    }

    std::vector<G4AttValue>* CreateAttValues() const override
    {
      auto* values = new std::vector<G4AttValue>;
      values->emplace_back("Pos", G4BestUnit(position_, "Length"), "");
      values->emplace_back("T", G4BestUnit(time_, "Time"), "");
      values->emplace_back("KE", G4BestUnit(kinetic_energy_, "Energy"), "");
      values->emplace_back("EDep", G4BestUnit(edep_, "Energy"), "");
      values->emplace_back("Proc", process_, "");
      values->emplace_back("Vol", volume_, "");
      return values;
    }

  private:
    G4ThreeVector position_;
    G4double time_;
    G4double kinetic_energy_;
    G4double edep_;
    G4String process_;
    G4String volume_;
};

// one replayed track, with the attributes of G4Trajectory so the standard
// vis models and filters (drawByParticleID, attributeFilter on IKE, ...) work
class RecordedTrajectory : public G4VTrajectory
{
  public:
    RecordedTrajectory(const RecordedTracks& tracks, std::size_t i)
        : track_id_(tracks.track_id[i]),
          parent_id_(tracks.parent_id[i]),
          pdg_(tracks.pdg[i]),
          particle_(tracks.particle[i]),
          charge_(tracks.charge[i]),
          creator_(tracks.creator_process[i])
    {
      const std::size_t first = tracks.track_point_start[i];
      const std::size_t end = i + 1 < tracks.track_point_start.size()
                                ? tracks.track_point_start[i + 1]
                                : tracks.x_mm.size();
      for (std::size_t p = first; p < end; ++p) {
        points_.push_back(std::make_unique<RecordedPoint>(tracks, p));
      }
      initial_kinetic_energy_ = tracks.ke_mev[first] * MeV;
    }

    G4int GetTrackID() const override { return track_id_; }
    G4int GetParentID() const override { return parent_id_; }
    G4String GetParticleName() const override { return particle_; }
    G4double GetCharge() const override { return charge_; }
    G4int GetPDGEncoding() const override { return pdg_; }
    G4ThreeVector GetInitialMomentum() const override { return {}; }
    G4int GetPointEntries() const override { return static_cast<G4int>(points_.size()); }
    G4VTrajectoryPoint* GetPoint(G4int i) const override
    {
      return points_[static_cast<std::size_t>(i)].get();
    }
    void AppendStep(const G4Step*) override {}
    void MergeTrajectory(G4VTrajectory*) override {}

    const std::map<G4String, G4AttDef>* GetAttDefs() const override
    {
      G4bool is_new;
      auto* store = G4AttDefStore::GetInstance("RecordedTrajectory", is_new);
      if (is_new) {
        (*store)["ID"] = G4AttDef("ID", "Track ID", "Physics", "", "G4int");
        (*store)["PID"] = G4AttDef("PID", "Parent ID", "Physics", "", "G4int");
        (*store)["PN"] = G4AttDef("PN", "Particle Name", "Physics", "", "G4String");
        (*store)["Ch"] = G4AttDef("Ch", "Charge", "Physics", "e+", "G4double");
        (*store)["PDG"] = G4AttDef("PDG", "PDG Encoding", "Physics", "", "G4int");
        (*store)["IKE"] =
          G4AttDef("IKE", "Initial kinetic energy", "Physics", "G4BestUnit", "G4double");
        (*store)["CPN"] = G4AttDef("CPN", "Creator Process Name", "Physics", "", "G4String");
        (*store)["NTP"] = G4AttDef("NTP", "No. of points", "Physics", "", "G4int");
      }
      return store;
    }

    std::vector<G4AttValue>* CreateAttValues() const override
    {
      auto* values = new std::vector<G4AttValue>;
      values->emplace_back("ID", G4UIcommand::ConvertToString(track_id_), "");
      values->emplace_back("PID", G4UIcommand::ConvertToString(parent_id_), "");
      values->emplace_back("PN", particle_, "");
      values->emplace_back("Ch", G4UIcommand::ConvertToString(charge_), "");
      values->emplace_back("PDG", G4UIcommand::ConvertToString(pdg_), "");
      values->emplace_back("IKE", G4BestUnit(initial_kinetic_energy_, "Energy"), "");
      values->emplace_back("CPN", creator_, "");
      values->emplace_back("NTP", G4UIcommand::ConvertToString(GetPointEntries()), "");
      return values;
    }

  private:
    G4int track_id_;
    G4int parent_id_;
    G4int pdg_;
    G4String particle_;
    G4double charge_;
    G4String creator_;
    G4double initial_kinetic_energy_ = 0.0;
    std::vector<std::unique_ptr<RecordedPoint>> points_;
};

// transport only: the viewer event exists to carry the recorded trajectories
class ViewerPhysicsList : public G4VUserPhysicsList
{
  public:
    void ConstructParticle() override
    {
      G4Geantino::Definition();
      G4Gamma::Definition();
      G4Electron::Definition();
      G4Positron::Definition();
      G4Proton::Definition();
    }
    void ConstructProcess() override { AddTransportation(); }
    void SetCuts() override { SetCutsWithDefault(); }
};

class GeantinoGenerator : public G4VUserPrimaryGeneratorAction
{
  public:
    GeantinoGenerator() : gun_(1)
    {
      gun_.SetParticleDefinition(G4Geantino::Definition());
      gun_.SetParticleEnergy(1.0 * MeV);
    }
    void GeneratePrimaries(G4Event* event) override { gun_.GeneratePrimaryVertex(event); }

  private:
    G4ParticleGun gun_;
};

class InjectTracksAction : public G4UserEventAction
{
  public:
    explicit InjectTracksAction(const RecordedTracks& tracks) : tracks_(tracks) {}

    void EndOfEventAction(const G4Event* event) override
    {
      auto* container = event->GetTrajectoryContainer();
      if (container == nullptr) {
        container = new G4TrajectoryContainer;
        const_cast<G4Event*>(event)->SetTrajectoryContainer(container);
      }
      for (std::size_t i = 0; i < tracks_.track_id.size(); ++i) {
        container->insert(new RecordedTrajectory(tracks_, i));
      }
    }

  private:
    const RecordedTracks& tracks_;
};

class ViewerActions : public G4VUserActionInitialization
{
  public:
    explicit ViewerActions(const RecordedTracks& tracks) : tracks_(tracks) {}
    void Build() const override
    {
      SetUserAction(new GeantinoGenerator);
      SetUserAction(new InjectTracksAction(tracks_));
    }

  private:
    const RecordedTracks& tracks_;
};

}  // namespace

void ViewEvent(
  const SessionConfig& config,
  const RecordedTracks& tracks,
  const std::vector<std::string>& commands)
{
  if (tracks.track_point_start.size() != tracks.track_id.size()) {
    throw std::runtime_error("Recorded tracks need one point start per track.");
  }

  static char program[] = "geant4_view";
  static char* argv[] = {program, nullptr};
  auto ui_executive = std::make_unique<G4UIExecutive>(1, argv, "qt");

  std::unique_ptr<G4RunManager> run_manager(
    G4RunManagerFactory::CreateRunManager(G4RunManagerType::SerialOnly));
  run_manager->SetUserInitialization(new DetectorConstruction(
    config.world_size_mm, config.detector_size_mm, config.envelope_material,
    config.device_components, config.em_production_cut_mm));
  run_manager->SetUserInitialization(new ViewerPhysicsList);
  run_manager->SetUserInitialization(new ViewerActions(tracks));
  run_manager->Initialize();

  auto vis_manager = std::make_unique<G4VisExecutive>("quiet");
  vis_manager->Initialize();

  auto* ui = G4UImanager::GetUIpointer();
  for (const auto& command : commands) {
    ui->ApplyCommand(command);
  }
  ui_executive->SessionStart();

  vis_manager.reset();
  run_manager.reset();
}

}  // namespace g4bridge
