#include "detector_construction.hpp"

#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "G4Box.hh"
#include "G4Element.hh"
#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4PVPlacement.hh"
#include "G4ProductionCuts.hh"
#include "G4Region.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"

namespace g4bridge
{

namespace
{

G4Material* BuildLiCoO2()
{
  const G4String name = "LiCoO2";
  if (auto* existing = G4Material::GetMaterial(name, false)) {
    return existing;
  }

  auto* nist = G4NistManager::Instance();
  auto* lithium = nist->FindOrBuildElement("Li");
  auto* cobalt = nist->FindOrBuildElement("Co");
  auto* oxygen = nist->FindOrBuildElement("O");
  if (!lithium || !cobalt || !oxygen) {
    throw std::runtime_error("Failed to build LiCoO2 elements.");
  }

  auto* material = new G4Material(name, 5.05 * g / cm3, 3);
  material->AddElement(lithium, 1);
  material->AddElement(cobalt, 1);
  material->AddElement(oxygen, 2);
  return material;
}

// Cured bisphenol-A epoxy resin (DGEBA, C21H24O4), the binder of mold
// compounds and FR-4.
G4Material* BuildEpoxy()
{
  const G4String name = "Epoxy";
  if (auto* existing = G4Material::GetMaterial(name, false)) {
    return existing;
  }

  auto* nist = G4NistManager::Instance();
  auto* material = new G4Material(name, 1.2 * g / cm3, 3);
  material->AddElement(nist->FindOrBuildElement("C"), 21);
  material->AddElement(nist->FindOrBuildElement("H"), 24);
  material->AddElement(nist->FindOrBuildElement("O"), 4);
  return material;
}

// Epoxy filled with silica, by mass: MoldCompound is a typical IC package
// compound (80% fused silica), FR4 a glass-epoxy circuit board (52.8% glass).
G4Material* BuildSilicaEpoxy(const G4String& name, G4double density, G4double silica_fraction)
{
  if (auto* existing = G4Material::GetMaterial(name, false)) {
    return existing;
  }

  auto* material = new G4Material(name, density, 2);
  material->AddMaterial(
    G4NistManager::Instance()->FindOrBuildMaterial("G4_SILICON_DIOXIDE"), silica_fraction);
  material->AddMaterial(BuildEpoxy(), 1.0 - silica_fraction);
  return material;
}

G4Material* ResolveMaterial(const std::string& name)
{
  if (name == "LiCoO2") {
    return BuildLiCoO2();
  }
  if (name == "MoldCompound") {
    return BuildSilicaEpoxy(name, 1.95 * g / cm3, 0.80);
  }
  if (name == "FR4") {
    return BuildSilicaEpoxy(name, 1.86 * g / cm3, 0.528);
  }

  auto* material = G4NistManager::Instance()->FindOrBuildMaterial(name);
  if (!material) {
    std::ostringstream msg;
    msg << "Failed to resolve Geant4 material '" << name << "'.";
    throw std::runtime_error(msg.str());
  }
  return material;
}

G4double HalfLength(double size_mm)
{
  return 0.5 * size_mm * mm;
}

G4ThreeVector PlacementOffset(
  const DeviceComponent& component,
  const DeviceComponent* parent)
{
  if (!parent) {
    return G4ThreeVector(
      component.center_mm[0] * mm,
      component.center_mm[1] * mm,
      component.center_mm[2] * mm);
  }
  return G4ThreeVector(
    (component.center_mm[0] - parent->center_mm[0]) * mm,
    (component.center_mm[1] - parent->center_mm[1]) * mm,
    (component.center_mm[2] - parent->center_mm[2]) * mm);
}

}  // namespace

DetectorConstruction::DetectorConstruction(
  const std::array<double, 3>& world_size_mm,
  const std::array<double, 3>& detector_size_mm,
  const std::string& envelope_material,
  std::vector<DeviceComponent> components,
  double em_production_cut_mm)
    : world_size_mm_(world_size_mm),
      detector_size_mm_(detector_size_mm),
      envelope_material_(envelope_material),
      components_(std::move(components)),
      em_production_cut_mm_(em_production_cut_mm)
{
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  constexpr G4bool check_overlaps = true;

  auto* world_mat = ResolveMaterial("G4_Galactic");
  auto* envelope_mat = ResolveMaterial(envelope_material_);

  auto* solid_world = new G4Box(
    "World",
    HalfLength(world_size_mm_[0]),
    HalfLength(world_size_mm_[1]),
    HalfLength(world_size_mm_[2]));
  auto* logic_world = new G4LogicalVolume(solid_world, world_mat, "World");

  auto* phys_world = new G4PVPlacement(
    nullptr, G4ThreeVector(), logic_world, "World", nullptr, false, 0, check_overlaps);

  auto* solid_envelope = new G4Box(
    "DetectorEnvelope",
    HalfLength(detector_size_mm_[0]),
    HalfLength(detector_size_mm_[1]),
    HalfLength(detector_size_mm_[2]));
  auto* logic_envelope =
    new G4LogicalVolume(solid_envelope, envelope_mat, "DetectorEnvelope");

  new G4PVPlacement(
    nullptr,
    G4ThreeVector(),
    logic_envelope,
    "DetectorEnvelope",
    logic_world,
    false,
    0,
    check_overlaps);

  scoring_volumes_.clear();
  scoring_names_.clear();
  electronics_region_ = nullptr;
  std::unordered_map<std::string, std::size_t> component_indices;
  component_indices.reserve(components_.size());
  for (std::size_t i = 0; i < components_.size(); ++i) {
    component_indices.emplace(components_[i].name, i);
  }

  std::vector<G4LogicalVolume*> logic_components(components_.size(), nullptr);
  for (std::size_t i = 0; i < components_.size(); ++i) {
    const auto& component = components_[i];
    auto* solid_component = new G4Box(
      component.name,
      HalfLength(component.size_mm[0]),
      HalfLength(component.size_mm[1]),
      HalfLength(component.size_mm[2]));
    auto* logic_component =
      new G4LogicalVolume(solid_component, ResolveMaterial(component.material), component.name);
    logic_components[i] = logic_component;
  }

  std::vector<bool> placed(components_.size(), false);
  auto place_component = [&](auto&& self, std::size_t i) -> void {
    if (placed[i]) {
      return;
    }
    const auto& component = components_[i];
    G4LogicalVolume* mother = logic_envelope;
    const DeviceComponent* parent = nullptr;
    if (!component.parent.empty()) {
      const auto parent_idx = component_indices.at(component.parent);
      self(self, parent_idx);
      mother = logic_components[parent_idx];
      parent = &components_[parent_idx];
    }

    new G4PVPlacement(
      nullptr,
      PlacementOffset(component, parent),
      logic_components[i],
      component.name,
      mother,
      false,
      static_cast<G4int>(i),
      check_overlaps);
    placed[i] = true;
  };

  for (std::size_t i = 0; i < components_.size(); ++i) {
    place_component(place_component, i);
  }
  for (std::size_t i = 0; i < components_.size(); ++i) {
    const auto& component = components_[i];
    if (component.score) {
      scoring_volumes_.push_back(logic_components[i]);
      scoring_names_.push_back(component.name);
    }
  }

  if (em_production_cut_mm_ > 0.0) {
    std::unordered_set<std::size_t> roots;
    for (std::size_t i = 0; i < components_.size(); ++i) {
      if (!components_[i].score) continue;
      const auto& parent = components_[i].parent;
      roots.insert(parent.empty() ? i : component_indices.at(parent));
    }
    if (!roots.empty()) {
      electronics_region_ = new G4Region("ElectronicsRegion");
      for (const auto root : roots) {
        bool nested = false;
        std::string parent = components_[root].parent;
        while (!parent.empty()) {
          const auto ancestor = component_indices.at(parent);
          if (roots.count(ancestor)) nested = true;
          parent = components_[ancestor].parent;
        }
        if (!nested) electronics_region_->AddRootLogicalVolume(logic_components[root]);
      }
      auto* cuts = new G4ProductionCuts();
      const G4double em_cut = em_production_cut_mm_ * mm;
      cuts->SetProductionCut(em_cut, "gamma");
      cuts->SetProductionCut(em_cut, "e-");
      cuts->SetProductionCut(em_cut, "e+");
      cuts->SetProductionCut(0.0, "proton");
      electronics_region_->SetProductionCuts(cuts);
    }
  }

  return phys_world;
}

}  // namespace g4bridge
