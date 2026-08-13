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

namespace g4bridge
{

BankPrimaryGeneratorAction::BankPrimaryGeneratorAction(
  const PrimaryBank& primary_bank,
  const SourceDistribution& source_distribution)
    : primary_bank_(primary_bank),
      source_distribution_(source_distribution),
      particle_gun_(std::make_unique<G4ParticleGun>(1))
{
}

void BankPrimaryGeneratorAction::GeneratePrimaries(G4Event* event)
{
  // choose active source for this event
  const std::size_t event_id = static_cast<std::size_t>(event->GetEventID());
  if (source_distribution_.Loaded()) {
    if (event_id >= source_distribution_.Size()) {
      std::ostringstream msg;
      msg << "Event " << event_id << " requested but source distribution size is "
          << source_distribution_.Size() << ".";
      throw std::runtime_error(msg.str());
    }
    GeneratePrimary(*event, source_distribution_.Sample());
    return;
  }

  if (event_id >= primary_bank_.Size()) {
    std::ostringstream msg;
    msg << "Event " << event_id << " requested but primary bank size is " << primary_bank_.Size()
        << ".";
    throw std::runtime_error(msg.str());
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
