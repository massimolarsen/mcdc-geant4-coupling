#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import os
import pathlib
import subprocess
import sys

import h5py
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


def _to_geant4_bank(handoff_particles: np.ndarray) -> np.ndarray:
    neutron_mask = handoff_particles["particle_type"] == 0
    if np.any(~neutron_mask):
        raise RuntimeError(
            "Handoff bank contains non-neutron particles; current Geant4 bridge supports neutrons only."
        )
    neutron_rows = handoff_particles[neutron_mask]

    bank = np.empty((len(neutron_rows), 10), dtype=np.float64)
    if len(neutron_rows) == 0:
        return bank

    bank[:, 0] = 2112.0
    bank[:, 1] = neutron_rows["x"] * 10.0
    bank[:, 2] = neutron_rows["y"] * 10.0
    bank[:, 3] = neutron_rows["z"] * 10.0
    bank[:, 4] = neutron_rows["ux"]
    bank[:, 5] = neutron_rows["uy"]
    bank[:, 6] = neutron_rows["uz"]
    bank[:, 7] = neutron_rows["E"] * 1.0e-6
    bank[:, 8] = neutron_rows["w"]
    bank[:, 9] = neutron_rows["t"] * 1.0e9
    return bank


def main() -> int:
    this_file = pathlib.Path(__file__).resolve()
    couple_repo = this_file.parents[2]
    repos_root = couple_repo.parent
    mcdc_example_dir = repos_root / "MCDC" / "examples" / "cube_in_cube"
    mcdc_input = mcdc_example_dir / "input_cell_current.py"
    output_path = mcdc_example_dir / "cube_in_cube_cell_current.h5"
    mcdc_lib = repos_root / "MCDC" / "hdf5lib"

    mcdc_env = dict(**os.environ)
    mcdc_env.setdefault("MCDC_LIB", str(mcdc_lib))
    if not pathlib.Path(mcdc_env["MCDC_LIB"]).exists():
        raise RuntimeError(
            f"MCDC_LIB path does not exist: {mcdc_env['MCDC_LIB']}. Set MCDC_LIB to your CE library directory."
        )

    subprocess.run(
        [sys.executable, str(mcdc_input), "--mode=numba"],
        cwd=mcdc_example_dir,
        check=True,
        env=mcdc_env,
    )

    with h5py.File(output_path, "r") as f:
        handoff_bank_size = int(f["handoff/particles_size"][()])
        handoff_particles = np.asarray(f["handoff/particles"][:handoff_bank_size])

    geant4_bank = _to_geant4_bank(handoff_particles)

    g4 = _load_local_geant4_bridge(couple_repo / "build")
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

    session = g4.Session(config)
    session.initialize()
    session.load_primaries(geant4_bank)
    session.beam_on()
    results = session.get_results()
    session.close()

    print(
        {
            "handoff_bank_size": handoff_bank_size,
            "loaded_primaries": int(results.loaded_primaries),
            "events_run": int(results.last_events_run),
            "match": int(results.loaded_primaries) == handoff_bank_size,
        }
    )

    if int(results.loaded_primaries) != handoff_bank_size:
        raise RuntimeError(
            "loaded_primaries does not match handoff bank size."
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
