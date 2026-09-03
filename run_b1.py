#!/usr/bin/env python3
from __future__ import annotations

import importlib
import pathlib
import sys

import numpy as np


def main() -> int:
    repo_dir = pathlib.Path(__file__).resolve().parent
    build_dir = repo_dir / "build"

    # Import the in-tree extension module built by CMake directly from build/.
    sys.path.insert(0, str(build_dir))
    g4 = importlib.import_module("geant4_bridge")

    config = g4.SessionConfig()
    config.world_size_mm = [100.0, 100.0, 100.0]
    config.detector_size_mm = [10.0, 10.0, 10.0]
    config.detector_material = "G4_Si"
    config.envelope_material = "G4_Galactic"
    config.physics_list = "QGSP_BIC"
    detector = g4.DeviceComponent()
    detector.name = "detector"
    detector.material = "G4_Si"
    detector.center_mm = [0.0, 0.0, 0.0]
    detector.size_mm = [10.0, 10.0, 10.0]
    detector.score = True
    config.device_components = [detector]

    bank = np.array(
        [
            # particle_id, x_mm, y_mm, z_mm, ux, uy, uz, E_MeV, weight, time_ns
            [2112, 0.0, 0.0, -40.0, 0.0, 0.0, 1.0, 14.0, 1.0, 0.0],
            [2112, 1.0, 0.0, -40.0, 0.0, 0.0, 1.0, 14.0, 2.0, 0.0],
        ],
        dtype=np.float64,
    )

    session = g4.Session(config)
    session.initialize()
    session.load_primaries(bank)
    session.beam_on()
    results = session.get_results()
    print(
        {
            "physics_list": results.physics_list,
            "loaded_primaries": results.loaded_primaries,
            "last_events_run": results.last_events_run,
            "last_total_edep_mev": results.last_total_edep_mev,
            "last_dose_gy": results.last_dose_gy,
        }
    )
    session.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
