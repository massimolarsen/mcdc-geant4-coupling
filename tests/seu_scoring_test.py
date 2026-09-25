"""Focused Geant4 bridge check for distribution-mode SEU scores."""

import pathlib
import sys
import tempfile

import numpy as np

sys.path.insert(0, sys.argv[1])
import geant4_bridge as bridge


def main():
    with tempfile.TemporaryDirectory() as diagnostic_dir:
        config = bridge.SessionConfig()
        config.world_size_mm = [100.0, 100.0, 100.0]
        config.detector_size_mm = [30.0, 30.0, 30.0]
        config.physics_list = "QGSP_BIC_HP"
        config.em_production_cut_mm = 0.001
        config.record_seu_events = True
        config.diagnostic_min_Eion_mev = 0.001
        config.diagnostic_dir = diagnostic_dir
        config.random_seed = 12345
        components = []
        for name, parent, size, score in (
            ("die", "", 20.0, False),
            ("sv_a", "die", 8.0, True),
            ("sv_b", "die", 8.0, True),
        ):
            component = bridge.DeviceComponent()
            component.name = name
            component.parent = parent
            component.material = "G4_Si"
            component.size_mm = [size, size, size]
            component.score = score
            component.center_mm = [-5.0, 0.0, 0.0] if name == "sv_a" else (
                [5.0, 0.0, 0.0] if name == "sv_b" else [0.0, 0.0, 0.0]
            )
            components.append(component)
        config.device_components = components

        session = bridge.Session(config)
        session.initialize()
        weights = np.zeros(6)
        weights[4] = 20.0
        session.load_source_distribution(
            np.array([[-8.0, 8.0], [-3.0, 3.0], [-3.0, 3.0]]),
            np.array([0.98, 1.0]),
            np.array([-np.pi, np.pi]),
            np.array([13.0, 15.0]),
            weights,
            1,
            1,
            100,
        )
        session.beam_on()
        result = session.get_results()
        session.close()

        assert result.last_events_run == 100
        assert result.em_production_cut_mm == 0.001
        assert result.proton_production_cut_mm == 0.0
        assert result.electronics_cut_materials == ["G4_Si"], result.electronics_cut_materials
        assert len(result.electronics_cut_energy_mev) == 4
        assert all(x > 0.0 for x in result.electronics_cut_energy_mev[:3])
        assert result.electronics_cut_energy_mev[3] == 0.0
        assert len(result.component_names) == 2
        assert sum(result.component_event_ionizing_count) == 200
        assert abs(sum(result.component_event_ionizing_sumw) - 40.0) < 1e-10
        assert abs(sum(result.component_event_ionizing_sumw2) - 8.0) < 1e-10
        for i in range(2):
            total = result.component_edep_mev[i]
            niel = result.component_niel_mev[i]
            ion = result.component_ionizing_mev[i]
            species = result.component_species_ionizing_mev[i * 7:(i + 1) * 7]
            assert abs(total - niel - ion) < 1e-8
            assert abs(sum(species) - ion) < 1e-8
            assert abs(result.component_primary_ionizing_mev[i] +
                       result.component_secondary_ionizing_mev[i] - ion) < 1e-8
        lines = pathlib.Path(diagnostic_dir, "worker_-1.tsv").read_text().splitlines()
        selected = [line for line in lines if line.startswith("M\t")]
        assert 0 < len(selected) < 100
        assert sum(line.startswith("E\t") for line in lines) == 2 * len(selected)
        assert any(line.startswith("I\t") for line in lines)
        assert any(line.startswith("I\t") and line.split("\t")[4] == "1"
                   for line in lines)
        assert {tuple(line.split("\t")[1:3]) for line in lines} == {
            tuple(line.split("\t")[1:3]) for line in selected
        }


if __name__ == "__main__":
    main()
