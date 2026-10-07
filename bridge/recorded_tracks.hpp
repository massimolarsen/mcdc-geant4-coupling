#ifndef G4BRIDGE_RECORDED_TRACKS_HPP
#define G4BRIDGE_RECORDED_TRACKS_HPP

#include <cstddef>
#include <string>
#include <vector>

namespace g4bridge
{

// Every track of one event with its step points, in tracking order. Track i
// owns points [track_point_start[i], track_point_start[i + 1]) (or to the end).
// Point 0 of a track is its vertex; later points end one step each, carrying
// that step's deposit, defining process, and the volume it was taken in.
struct RecordedTracks
{
  std::vector<int> track_id;
  std::vector<int> parent_id;
  std::vector<int> pdg;
  std::vector<std::string> particle;
  std::vector<double> charge;
  std::vector<std::string> creator_process;
  std::vector<std::size_t> track_point_start;

  std::vector<double> x_mm;
  std::vector<double> y_mm;
  std::vector<double> z_mm;
  std::vector<double> t_ns;
  std::vector<double> ke_mev;
  std::vector<double> edep_mev;
  std::vector<std::string> process;
  std::vector<std::string> volume;

  void Clear()
  {
    track_id.clear();
    parent_id.clear();
    pdg.clear();
    particle.clear();
    charge.clear();
    creator_process.clear();
    track_point_start.clear();
    x_mm.clear();
    y_mm.clear();
    z_mm.clear();
    t_ns.clear();
    ke_mev.clear();
    edep_mev.clear();
    process.clear();
    volume.clear();
  }

  void AddPoint(
    double x, double y, double z, double t, double ke, double edep, std::string process_name,
    std::string volume_name)
  {
    x_mm.push_back(x);
    y_mm.push_back(y);
    z_mm.push_back(z);
    t_ns.push_back(t);
    ke_mev.push_back(ke);
    edep_mev.push_back(edep);
    process.push_back(std::move(process_name));
    volume.push_back(std::move(volume_name));
  }
};

}  // namespace g4bridge

#endif
