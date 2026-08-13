#ifndef G4BRIDGE_DEVICE_GEOMETRY_HPP
#define G4BRIDGE_DEVICE_GEOMETRY_HPP

#include <array>
#include <string>

namespace g4bridge
{

struct DeviceComponent
{
  std::string name = "";
  std::string material = "";
  std::string parent = "";
  std::array<double, 3> center_mm{0.0, 0.0, 0.0};
  std::array<double, 3> size_mm{1.0, 1.0, 1.0};
  bool score = true;
};

}  // namespace g4bridge

#endif
