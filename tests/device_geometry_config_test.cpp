#include "session.hpp"

#include <iostream>
#include <stdexcept>
#include <string>

namespace
{

void Require(bool condition, const char* message)
{
  if (!condition) {
    throw std::runtime_error(message);
  }
}

g4bridge::SessionConfig BaseConfig()
{
  g4bridge::SessionConfig config;
  config.world_size_mm = {24.0, 24.0, 24.0};
  config.detector_size_mm = {20.0, 20.0, 20.0};
  config.envelope_material = "G4_Galactic";

  g4bridge::DeviceComponent component;
  component.name = "detector";
  component.material = "G4_Si";
  component.center_mm = {0.0, 0.0, 0.0};
  component.size_mm = {20.0, 20.0, 20.0};
  component.score = true;
  config.device_components = {component};
  return config;
}

void ExpectThrows(const g4bridge::SessionConfig& config, const std::string& match)
{
  try {
    g4bridge::Session session(config);
  }
  catch (const std::exception& exc) {
    Require(std::string(exc.what()).find(match) != std::string::npos, exc.what());
    return;
  }
  throw std::runtime_error("Expected Session constructor to reject invalid config.");
}

}  // namespace

int main()
{
  try {
    {
      auto config = BaseConfig();
      config.device_components.clear();
      ExpectThrows(config, "device_components");
    }
    {
      auto config = BaseConfig();
      config.device_components[0].size_mm = {0.0, 1.0, 1.0};
      ExpectThrows(config, "positive");
    }
    {
      auto config = BaseConfig();
      config.n_threads = 0;
      ExpectThrows(config, "n_threads");
    }
    {
      auto config = BaseConfig();
      config.device_components[0].center_mm = {1.0, 0.0, 0.0};
      ExpectThrows(config, "fit");
    }
    {
      auto config = BaseConfig();
      auto second = config.device_components[0];
      second.name = "overlap";
      second.size_mm = {2.0, 2.0, 2.0};
      config.device_components.push_back(second);
      ExpectThrows(config, "overlap");
    }
    {
      auto config = BaseConfig();
      auto child = config.device_components[0];
      child.name = "child";
      child.parent = "detector";
      child.size_mm = {2.0, 2.0, 2.0};
      config.device_components.push_back(child);
      g4bridge::Session session(config);
    }
    {
      auto config = BaseConfig();
      config.device_components[0].parent = "missing";
      ExpectThrows(config, "unknown parent");
    }
    {
      auto config = BaseConfig();
      config.device_components[0].parent = "detector";
      ExpectThrows(config, "parent itself");
    }
    {
      auto config = BaseConfig();
      auto child = config.device_components[0];
      child.name = "child";
      child.parent = "detector";
      config.device_components[0].parent = "child";
      config.device_components.push_back(child);
      ExpectThrows(config, "cycle");
    }
    {
      auto config = BaseConfig();
      auto child = config.device_components[0];
      child.name = "child";
      child.parent = "detector";
      child.size_mm = {4.0, 4.0, 4.0};
      child.center_mm = {9.0, 0.0, 0.0};
      config.device_components.push_back(child);
      ExpectThrows(config, "parent");
    }
    {
      auto config = BaseConfig();
      auto child_a = config.device_components[0];
      child_a.name = "child_a";
      child_a.parent = "detector";
      child_a.size_mm = {4.0, 4.0, 4.0};
      child_a.center_mm = {-1.0, 0.0, 0.0};
      auto child_b = child_a;
      child_b.name = "child_b";
      child_b.center_mm = {1.0, 0.0, 0.0};
      config.device_components.push_back(child_a);
      config.device_components.push_back(child_b);
      ExpectThrows(config, "overlap");
    }
    g4bridge::Session session(BaseConfig());
    return 0;
  }
  catch (const std::exception& exc) {
    std::cerr << exc.what() << '\n';
    return 1;
  }
}
