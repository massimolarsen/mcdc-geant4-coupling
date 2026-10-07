#ifndef G4BRIDGE_EVENT_VIEWER_HPP
#define G4BRIDGE_EVENT_VIEWER_HPP

#include <string>
#include <vector>

#include "recorded_tracks.hpp"
#include "session.hpp"

namespace g4bridge
{

// Open Geant4's Qt viewer on the geometry in config, showing the recorded
// tracks as the trajectories of one event. commands run before the UI session
// starts and must include the /run/beamOn 1 that draws the event.
void ViewEvent(
  const SessionConfig& config,
  const RecordedTracks& tracks,
  const std::vector<std::string>& commands);

}  // namespace g4bridge

#endif
