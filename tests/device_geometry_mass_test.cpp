#include "detector_construction.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4SystemOfUnits.hh"

namespace
{

void RequireClose(G4double actual, G4double expected, const char* name)
{
  const G4double scale = std::max(std::abs(expected), 1.0e-30);
  if (std::abs(actual - expected) / scale > 1.0e-12) {
    std::cerr << name << " actual=" << actual / kg
              << " kg expected=" << expected / kg << " kg\n";
    throw std::runtime_error("Mass check failed.");
  }
}

G4double BoxVolume(double x_mm, double y_mm, double z_mm)
{
  return (x_mm * mm) * (y_mm * mm) * (z_mm * mm);
}

}  // namespace

int main()
{
  try {
    g4bridge::DeviceComponent package;
    package.name = "package";
    package.material = "G4_BAKELITE";
    package.size_mm = {10.0, 10.0, 10.0};
    package.score = true;

    g4bridge::DeviceComponent die;
    die.name = "die";
    die.material = "G4_Si";
    die.parent = "package";
    die.size_mm = {4.0, 4.0, 4.0};
    die.score = true;

    g4bridge::DeviceComponent active;
    active.name = "active";
    active.material = "G4_Si";
    active.parent = "die";
    active.size_mm = {2.0, 2.0, 2.0};
    active.score = true;

    g4bridge::DetectorConstruction detector(
      {24.0, 24.0, 24.0},
      {20.0, 20.0, 20.0},
      "G4_Galactic",
      {package, die, active});
    detector.Construct();

    const auto& scoring_volumes = detector.GetScoringVolumes();
    if (scoring_volumes.size() != 3) {
      throw std::runtime_error("Unexpected scoring volume count.");
    }

    const G4double package_mass =
      scoring_volumes[0]->GetMaterial()->GetDensity() *
      (BoxVolume(10.0, 10.0, 10.0) - BoxVolume(4.0, 4.0, 4.0));
    const G4double die_mass =
      scoring_volumes[1]->GetMaterial()->GetDensity() *
      (BoxVolume(4.0, 4.0, 4.0) - BoxVolume(2.0, 2.0, 2.0));
    const G4double active_mass =
      scoring_volumes[2]->GetMaterial()->GetDensity() * BoxVolume(2.0, 2.0, 2.0);

    RequireClose(scoring_volumes[0]->GetMass(false, false), package_mass, "package");
    RequireClose(scoring_volumes[1]->GetMass(false, false), die_mass, "die");
    RequireClose(scoring_volumes[2]->GetMass(false, false), active_mass, "active");
    return 0;
  }
  catch (const std::exception& exc) {
    std::cerr << exc.what() << '\n';
    return 1;
  }
}
