//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file DetectorConstruction.cc
/// \brief Implementation of the B1::DetectorConstruction class

#include "DetectorConstruction.hh"

#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"

#include <sstream>
#include <stdexcept>

namespace B1
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction()
    : DetectorConstruction({100.0, 100.0, 100.0}, {10.0, 10.0, 10.0}, "G4_Si")
{
}

DetectorConstruction::DetectorConstruction(
  const std::array<double, 3>& world_size_mm,
  const std::array<double, 3>& detector_size_mm,
  const std::string& detector_material)
    : world_size_mm_(world_size_mm),
      detector_size_mm_(detector_size_mm),
      detector_material_(detector_material)
{
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();

  // Option to switch on/off checking of volumes overlaps
  G4bool checkOverlaps = true;

  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");

  G4Material* detector_mat = nist->FindOrBuildMaterial(detector_material_);
  if (!detector_mat) {
    std::ostringstream msg;
    msg << "Failed to resolve detector material '" << detector_material_ << "'.";
    throw std::runtime_error(msg.str());
  }

  const G4double world_hx = 0.5 * world_size_mm_[0] * mm;
  const G4double world_hy = 0.5 * world_size_mm_[1] * mm;
  const G4double world_hz = 0.5 * world_size_mm_[2] * mm;

  const G4double detector_hx = 0.5 * detector_size_mm_[0] * mm;
  const G4double detector_hy = 0.5 * detector_size_mm_[1] * mm;
  const G4double detector_hz = 0.5 * detector_size_mm_[2] * mm;

  auto solidWorld = new G4Box("World", world_hx, world_hy, world_hz);
  auto logicWorld = new G4LogicalVolume(solidWorld, world_mat, "World");

  auto physWorld = new G4PVPlacement(
    nullptr, G4ThreeVector(), logicWorld, "World", nullptr, false, 0, checkOverlaps);

  auto solidDetector = new G4Box("Detector", detector_hx, detector_hy, detector_hz);
  auto logicDetector = new G4LogicalVolume(solidDetector, detector_mat, "Detector");

  new G4PVPlacement(
    nullptr, G4ThreeVector(), logicDetector, "Detector", logicWorld, false, 0, checkOverlaps);

  fScoringVolume = logicDetector;

  //
  // always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}  // namespace B1
