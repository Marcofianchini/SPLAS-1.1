import tempfile
import unittest
from pathlib import Path

import numpy as np

from splas.analysis import analyze_results
from splas.dataset import write_dataset
from splas.figures import make_run_figures


class FigurePipelineTests(unittest.TestCase):
    def test_synthetic_dataset_analysis_and_figures(self):
        time = np.arange(4.0)
        depth = np.array([2.5, 7.5, 12.5])
        phyto_size = np.geomspace(1.0, 200.0, 4)
        zoo_size = np.geomspace(2.0, 1000.0, 4)
        p = np.full((4, 3, 4), 0.01)
        z = np.full((4, 3, 4), 0.001)
        quota = np.full((4, 3, 4), 1e-12)
        light = 100.0 * np.exp(-0.1 * depth)[None, :] * np.ones((4, 1))
        coords = {
            "time": time,
            "depth": depth,
            "phyto_size": phyto_size,
            "zoo_size": zoo_size,
            "nutrient": np.arange(3.0),
        }
        data = {
            "P": p,
            "Z": z,
            "quota_P": quota,
            "cell_number_P": p / quota,
            "nutrient_pool": np.ones((4, 3, 3)),
            "detritus": np.zeros((4, 3)),
            "surface_irradiance": np.full(4, 100.0),
            "irradiance": light,
            "nutrient_limitation": np.full((4, 3, 4), 0.2),
            "grazing_limitation": np.full((4, 3, 4), 0.1),
            "light_limitation": np.full((4, 3, 4), 0.8),
            "phyto_growth_rate": np.full((4, 3, 4), 0.5),
            "grazing_loss": np.full((4, 3, 4), 1e-4),
            "net_phyto_growth_rate": np.full((4, 3, 4), 0.39),
            "total_growth_reduction": np.full((4, 3, 4), 0.61),
            "phyto_biomass": p.sum(axis=-1),
            "zoo_biomass": z.sum(axis=-1),
            "shannon_entropy_P": np.full((4, 3), np.log(4.0)),
            "shannon_evenness_P": np.ones((4, 3)),
            "diversity_valid_P": np.ones((4, 3)),
            "carbon_per_cell_P": np.ones(4),
            "carbon_biomass_P": np.ones((4, 3, 4)),
            "export_phyto": np.full(4, 1.0),
            "export_detritus": np.full(4, 2.0),
            "export_total": np.full(4, 3.0),
            "mass_balance_error": np.zeros(4),
            "nonnegative_correction": np.zeros(4),
        }
        config = {
            "name": "synthetic",
            "analysis": {"transient_days": 0.0},
            "forcing": {"multiplier": 1.0},
            "system": {"mode": "open"},
        }
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            run_dir = root / "run"
            run_dir.mkdir()
            write_dataset(
                run_dir / "results.nc",
                coords,
                data,
                config,
                {"config_hash": "synthetic"},
            )
            summary = analyze_results(run_dir, root / "tables")
            figures = make_run_figures(run_dir, root / "figures")
            self.assertEqual(summary["run_name"], "synthetic")
            self.assertEqual(len(figures), 6)
            self.assertTrue(all(path.stat().st_size > 1000 for path in figures))
