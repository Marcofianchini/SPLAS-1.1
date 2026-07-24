import tempfile
import unittest
from pathlib import Path

import numpy as np

from splas.config import ConfigError, config_hash, load_config, validate_config
from splas.dataset import read_dataset, write_dataset


ROOT = Path(__file__).resolve().parents[1]


class ConfigTests(unittest.TestCase):
    def test_all_versioned_configs_validate(self):
        configs = [
            path
            for path in (ROOT / "configs").rglob("*.json")
            if "manifests" not in path.parts
        ]
        self.assertGreaterEqual(len(configs), 30)
        for path in configs:
            with self.subTest(path=path):
                load_config(path)

    def test_closed_mode_rejects_transport(self):
        config = load_config(ROOT / "configs" / "reference" / "reference_closed_constant_default.json")
        config["system"]["diffusion_m2_s"] = 1e-4
        with self.assertRaises(ConfigError):
            validate_config(config)

    def test_hash_ignores_source_path(self):
        config = load_config(ROOT / "configs" / "reference" / "smoke.json")
        first = config_hash(config)
        config["_source"] = "/different/path.json"
        self.assertEqual(first, config_hash(config))


class DatasetTests(unittest.TestCase):
    def test_round_trip_self_describing_dataset(self):
        coords = {
            "time": np.array([0.0, 1.0]),
            "depth": np.array([2.5]),
            "phyto_size": np.array([1.0]),
            "zoo_size": np.array([2.0]),
            "nutrient": np.array([0.0, 1.0, 2.0]),
        }
        data = {
            "P": np.ones((2, 1, 1)),
            "Z": np.ones((2, 1, 1)),
            "quota_P": np.ones((2, 1, 1)),
            "cell_number_P": np.ones((2, 1, 1)),
            "nutrient_pool": np.ones((2, 1, 3)),
            "detritus": np.ones((2, 1)),
            "surface_irradiance": np.ones(2),
            "irradiance": np.ones((2, 1)),
            "nutrient_limitation": np.ones((2, 1, 1)),
            "grazing_limitation": np.ones((2, 1, 1)),
            "light_limitation": np.ones((2, 1, 1)),
            "phyto_growth_rate": np.ones((2, 1, 1)),
            "grazing_loss": np.ones((2, 1, 1)),
            "net_phyto_growth_rate": np.ones((2, 1, 1)),
            "total_growth_reduction": np.ones((2, 1, 1)),
            "phyto_biomass": np.ones((2, 1)),
            "zoo_biomass": np.ones((2, 1)),
            "shannon_entropy_P": np.ones((2, 1)),
            "shannon_evenness_P": np.ones((2, 1)),
            "diversity_valid_P": np.ones((2, 1)),
            "carbon_per_cell_P": np.ones(1),
            "carbon_biomass_P": np.ones((2, 1, 1)),
            "export_phyto": np.ones(2),
            "export_detritus": np.ones(2),
            "export_total": np.ones(2),
            "mass_balance_error": np.zeros(2),
            "nonnegative_correction": np.zeros(2),
        }
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / "results.nc"
            write_dataset(target, coords, data, {"name": "test"}, {"config_hash": "abc"})
            arrays, attrs = read_dataset(target)
        np.testing.assert_array_equal(arrays["P"], data["P"])
        self.assertEqual(attrs["config"]["name"], "test")
        self.assertEqual(attrs["metadata"]["config_hash"], "abc")
