"""Run separately with 1 and 2 threads to check serial and worker merging."""

import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(sys.argv[2]).resolve()))
import geant4_bridge as g4

config = g4.SessionConfig()
config.n_threads = int(sys.argv[1])
config.world_size_mm = [100.0, 100.0, 100.0]
config.detector_size_mm = [50.0, 50.0, 50.0]
config.detector_material = "G4_Si"
config.envelope_material = "G4_Galactic"
config.physics_list = "QGSP_BIC"
detector = g4.DeviceComponent()
detector.name = "test_silicon"
detector.material = "G4_Si"
detector.center_mm = [0.0, 0.0, 0.0]
detector.size_mm = [50.0, 50.0, 50.0]
detector.score = True
config.device_components = [detector]

# Charged particles starting inside silicon reliably exercise nonzero scoring.
bank = np.tile([11, 0, 0, 0, 0, 0, 1, 1.0, 1.0, 0], (40, 1)).astype(float)
session = g4.Session(config)
session.initialize()
for _ in range(2):
    session.load_primaries(bank)
    session.beam_on()
    result = session.get_results()
    assert result.last_events_run == len(bank)
    assert sum(result.edep_spectrum_counts) == len(bank)
    assert result.last_total_edep_mev > 0
    assert np.isclose(sum(result.edep_spectrum_edep_mev), result.last_total_edep_mev)
session.close()
print(f"Scoring and reset passed with {config.n_threads} thread(s)")
