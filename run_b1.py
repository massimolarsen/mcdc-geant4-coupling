#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import pathlib

import numpy as np


def _load_local_geant4_bridge(build_dir: pathlib.Path):
    extension_candidates = sorted(build_dir.glob("geant4_bridge*.so"))
    if not extension_candidates:
        raise RuntimeError(
            f"Could not find built geant4_bridge extension under: {build_dir}"
        )

    module_path = extension_candidates[0]
    spec = importlib.util.spec_from_file_location("geant4_bridge", module_path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Failed to load module spec from: {module_path}")

    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main() -> int:
    repo_dir = pathlib.Path(__file__).resolve().parent
    build_dir = repo_dir / "build"

    # Import the in-tree extension module built by CMake directly from build/.
    g4 = _load_local_geant4_bridge(build_dir)

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
            "initialized": results.initialized,
            "has_source": results.has_source,
            "physics_list": results.physics_list,
            "loaded_primaries": results.loaded_primaries,
            "last_events_run": results.last_events_run,
            "last_total_edep_mev": results.last_total_edep_mev,
            "last_dose_gy": results.last_dose_gy,
            "status": results.status,
        }
    )
    session.close()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
