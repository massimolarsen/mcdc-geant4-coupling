#ifndef G4BRIDGE_SEU_SCORING_HPP
#define G4BRIDGE_SEU_SCORING_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <string>
#include <vector>

namespace g4bridge
{

constexpr std::size_t kSeUSpeciesCount = 7;
constexpr std::size_t kSeUBins = 150;
constexpr std::size_t kSeUStride = kSeUBins + 3;
constexpr double kSeUMinMeV = 0.001;
constexpr double kSeUMaxMeV = 50.0;

inline std::vector<std::string> SeUSpeciesNames()
{
  return {"electron_positron", "proton", "neutron", "gamma", "alpha", "ion_recoil", "other"};
}

inline std::size_t SeUSpeciesForPdg(int pdg)
{
  const int code = std::abs(pdg);
  if (code == 11) return 0;
  if (code == 2212) return 1;
  if (code == 2112) return 2;
  if (code == 22) return 3;
  if (code == 1000020040) return 4;
  if (code >= 1000000000) return 5;
  return 6;
}

inline std::vector<double> SeUEdgesMeV()
{
  std::vector<double> edges(kSeUBins + 1);
  for (std::size_t i = 0; i <= kSeUBins; ++i) {
    edges[i] = kSeUMinMeV *
               std::pow(kSeUMaxMeV / kSeUMinMeV, double(i) / double(kSeUBins));
  }
  return edges;
}

// Slots: zero, underflow, logarithmic bins, overflow.
inline std::size_t SeUBin(double energy_mev, const std::vector<double>& edges)
{
  if (!(energy_mev > 0.0)) return 0;
  if (energy_mev < kSeUMinMeV) return 1;
  if (energy_mev >= kSeUMaxMeV) return kSeUStride - 1;
  return 2 + std::size_t(std::upper_bound(edges.begin(), edges.end(), energy_mev) -
                         edges.begin() - 1);
}

}  // namespace g4bridge

#endif
