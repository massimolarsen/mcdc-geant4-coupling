#include "bank_primary_generator_action.hpp"

#include <sstream>
#include <stdexcept>

#include "primary_bank.hpp"
#include "source_distribution.hpp"

#include "G4Event.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"

namespace g4bridge
{

BankPrimaryGeneratorAction::BankPrimaryGeneratorAction(
  const PrimaryBank& primary_bank,
  const SourceDistribution& source_distribution,
  int event_id_offset,
  std::vector<unsigned long> replay_state,
  bool capture_rng_state)
    : primary_bank_(primary_bank),
      source_distribution_(source_distribution),
      particle_gun_(std::make_unique<G4ParticleGun>(1)),
      event_id_offset_(event_id_offset),
      replay_state_(std::move(replay_state)),
      capture_rng_state_(capture_rng_state)
{
}

void BankPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  // Geant4 has just seeded this event; nothing has drawn from it yet
  if (!replay_state_.empty() && !G4Random::getTheEngine()->get(replay_state_)) {
    throw std::runtime_error("Replay RNG state does not match the Geant4 random engine.");
  }
  if (capture_rng_state_) {
    rng_state_ = G4Random::getTheEngine()->put();
  }

  // a replayed event keeps the ID it had in its original run
  if (event_id_offset_ != 0) {
    event->SetEventID(event->GetEventID() + event_id_offset_);
  }

  // choose active source for this event
  const std::size_t event_id = static_cast<std::size_t>(event->GetEventID());
  if (source_distribution_.Loaded()) {
    GeneratePrimary(*event, source_distribution_.Sample(event_id));
    return;
  }

  GeneratePrimary(*event, primary_bank_.At(event_id));
}

void BankPrimaryGeneratorAction::GeneratePrimary(
  G4Event& event,
  const Primary& p) const
{
  // look up particle type
  G4ParticleDefinition* particle =
    G4ParticleTable::GetParticleTable()->FindParticle(p.particle_id);
  if (!particle) {
    std::ostringstream msg;
    msg << "Failed to resolve particle ID " << p.particle_id << " in Geant4 particle table.";
    throw std::runtime_error(msg.str());
  }

  // load primary values into the Geant4 particle gun
  particle_gun_->SetParticleDefinition(particle);
  particle_gun_->SetParticlePosition(G4ThreeVector(p.x_mm * mm, p.y_mm * mm, p.z_mm * mm));
  particle_gun_->SetParticleMomentumDirection(G4ThreeVector(p.ux, p.uy, p.uz));
  particle_gun_->SetParticleEnergy(p.energy_mev * MeV);
  particle_gun_->SetParticleTime(p.time_ns * ns);
  particle_gun_->GeneratePrimaryVertex(&event);

  // preserve MCDC particle weight on the primary only. Geant4 forms the track
  // weight from vertex weight * primary weight, so setting both squares it.
  G4PrimaryVertex* vertex = event.GetPrimaryVertex(0);
  vertex->GetPrimary()->SetWeight(p.weight);

}

}  // namespace g4bridge
