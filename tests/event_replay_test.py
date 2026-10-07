"""Replaying one event, by seeds or by saved RNG state, reproduces it exactly."""

import json
import pathlib
import subprocess
import sys
import tempfile

import numpy as np

BUILD_DIR = sys.argv[1]
N_EVENTS = 30
SEED = 12345


def run(mode, diagnostic_dir, event=0, state=None):
    sys.path.insert(0, BUILD_DIR)
    import geant4_bridge as bridge

    config = bridge.SessionConfig()
    config.world_size_mm = [100.0, 100.0, 100.0]
    config.detector_size_mm = [30.0, 30.0, 30.0]
    config.physics_list = "QGSP_BIC_HP"
    config.em_production_cut_mm = 0.001
    config.random_seed = SEED
    config.record_seu_events = True
    config.diagnostic_min_Eion_mev = 0.0
    config.diagnostic_dir = diagnostic_dir
    config.rng_state_min_Eion_mev = 1.0e-9
    components = []
    for name, parent, size, center, score in (
        ("die", "", 20.0, [0.0, 0.0, 0.0], False),
        ("sv_a", "die", 8.0, [-5.0, 0.0, 0.0], True),
        ("sv_b", "die", 8.0, [5.0, 0.0, 0.0], True),
    ):
        component = bridge.DeviceComponent()
        component.name = name
        component.parent = parent
        component.material = "G4_Si"
        component.size_mm = [size, size, size]
        component.center_mm = center
        component.score = score
        components.append(component)
    config.device_components = components
    if mode == "original":
        config.n_threads = 2
    else:
        config.n_threads = 1
        config.replay_mode = mode
        config.event_id_offset = event
        config.record_tracks = True
        if state is not None:
            config.replay_state = state

    session = bridge.Session(config)
    session.initialize()
    weights = np.zeros(6)
    weights[4] = 20.0
    session.load_source_distribution(
        np.array([[-8.0, 8.0], [-3.0, 3.0], [-3.0, 3.0]]),
        np.array([0.98, 1.0]),
        np.array([-np.pi, np.pi]),
        np.array([50.0, 60.0]),
        weights,
        1,
        1,
        N_EVENTS,
        2212,
    )
    session.beam_on(-1 if mode == "original" else 1)
    tracks = session.get_results().tracks
    session.close()
    if mode != "original":
        # step deposits per volume, to check against the scored deposit
        edep = {}
        for volume, value in zip(tracks.volume, tracks.edep_mev):
            edep[volume] = edep.get(volume, 0.0) + value
        print(json.dumps({"n_tracks": len(tracks.track_id), "edep": edep}))


def rows(diagnostic_dir):
    lines = []
    for path in pathlib.Path(diagnostic_dir).glob("worker_*.tsv"):
        lines += path.read_text().splitlines()
    return lines


def event_rows(lines, event, codes="EMIN"):
    return sorted(
        line for line in lines if line[0] in codes and line.split("\t")[2] == str(event)
    )


def child(mode, event=0, state=None):
    with tempfile.TemporaryDirectory() as diagnostic_dir:
        args = [sys.executable, __file__, BUILD_DIR, mode, diagnostic_dir, str(event)]
        if state is not None:
            args.append(" ".join(state))
        out = subprocess.run(args, check=True, capture_output=True, text=True).stdout
        summary = [line for line in out.splitlines() if line.startswith("{")]
        return rows(diagnostic_dir), json.loads(summary[-1]) if summary else None


def main():
    original, _ = child("original")
    assert sum(line.startswith("M\t") for line in original) == N_EVENTS
    states = {line.split("\t")[2]: line.split("\t")[3].split() for line in original
              if line.startswith("R\t")}
    assert states, "no RNG states were written"

    # a late event, so the seed replay skips most of the master's draws
    event = max(int(e) for e in states)
    expected = event_rows(original, event)
    for mode, state in (("seeds", None), ("state", states[str(event)])):
        replayed, summary = child(mode, event, state)
        assert event_rows(replayed, event) == expected, f"{mode} replay of event {event} differs"
        assert event_rows(replayed, event, "R") == event_rows(original, event, "R")
        assert summary["n_tracks"] > 0
        for line in expected:
            fields = line.split("\t")
            if fields[0] == "E":
                volume = ("sv_a", "sv_b")[int(fields[3])]
                assert abs(summary["edep"].get(volume, 0.0) - float(fields[6])) < 1e-9


if __name__ == "__main__":
    if len(sys.argv) > 2:
        mode, diagnostic_dir, event = sys.argv[2:5]
        state = [int(x) for x in sys.argv[5].split()] if len(sys.argv) > 5 else None
        run(mode, diagnostic_dir, int(event), state)
    else:
        main()
