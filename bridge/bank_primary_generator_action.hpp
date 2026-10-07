#ifndef G4BRIDGE_BANK_PRIMARY_GENERATOR_ACTION_HPP
#define G4BRIDGE_BANK_PRIMARY_GENERATOR_ACTION_HPP

#include <cstddef>
#include <memory>
#include <vector>

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
      const SourceDistribution& source_distribution,
      int event_id_offset = 0,
      std::vector<unsigned long> replay_state = {},
      bool capture_rng_state = false);
    ~BankPrimaryGeneratorAction() override = default;

    void GeneratePrimaries(G4Event* event) override;

    // engine state at the start of the current event, if captured
    const std::vector<unsigned long>& EventRngState() const { return rng_state_; }

  private:
    // convert bridge primary data into a Geant4 vertex
    void GeneratePrimary(G4Event& event, const Primary& primary) const;

    const PrimaryBank& primary_bank_;
    const SourceDistribution& source_distribution_;
    std::unique_ptr<G4ParticleGun> particle_gun_;
    int event_id_offset_;
    std::vector<unsigned long> replay_state_;
    bool capture_rng_state_;
    std::vector<unsigned long> rng_state_;
};

}  // namespace g4bridge

#endif
