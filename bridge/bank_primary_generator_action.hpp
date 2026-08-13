#ifndef G4BRIDGE_BANK_PRIMARY_GENERATOR_ACTION_HPP
#define G4BRIDGE_BANK_PRIMARY_GENERATOR_ACTION_HPP

#include <cstddef>
#include <memory>

#include "G4VUserPrimaryGeneratorAction.hh"

class G4Event;
class G4ParticleGun;

namespace g4bridge
{

class PrimaryBank;
class SourceDistribution;
struct Primary;

class BankPrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
  public:
    // connect Geant4 primary generation to the active bridge source
    BankPrimaryGeneratorAction(
      const PrimaryBank& primary_bank,
      const SourceDistribution& source_distribution);
    ~BankPrimaryGeneratorAction() override = default;

    void GeneratePrimaries(G4Event* event) override;

  private:
    // convert bridge primary data into a Geant4 vertex
    void GeneratePrimary(G4Event& event, const Primary& primary) const;

    const PrimaryBank& primary_bank_;
    const SourceDistribution& source_distribution_;
    std::unique_ptr<G4ParticleGun> particle_gun_;
};

}  // namespace g4bridge

#endif
