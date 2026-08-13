#ifndef G4BRIDGE_DETECTOR_CONSTRUCTION_HPP
#define G4BRIDGE_DETECTOR_CONSTRUCTION_HPP

#include <array>
#include <string>
#include <vector>

#include "G4VUserDetectorConstruction.hh"

#include "device_geometry.hpp"

class G4LogicalVolume;
class G4VPhysicalVolume;

namespace g4bridge
{

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    DetectorConstruction(
      const std::array<double, 3>& world_size_mm,
      const std::array<double, 3>& detector_size_mm,
      const std::string& envelope_material,
      std::vector<DeviceComponent> components);
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;

    const std::vector<G4LogicalVolume*>& GetScoringVolumes() const
    {
      return scoring_volumes_;
    }

    const std::vector<std::string>& GetScoringNames() const
    {
      return scoring_names_;
    }

  private:
    std::array<double, 3> world_size_mm_;
    std::array<double, 3> detector_size_mm_;
    std::string envelope_material_;
    std::vector<DeviceComponent> components_;
    std::vector<G4LogicalVolume*> scoring_volumes_;
    std::vector<std::string> scoring_names_;
};

}  // namespace g4bridge

#endif
