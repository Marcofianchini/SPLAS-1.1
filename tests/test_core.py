import copy
import os
import unittest
from pathlib import Path

import numpy as np

from splas.backend import CoreUnavailable, require_core, run_core
from splas.config import load_config


ROOT = Path(__file__).resolve().parents[1]


def core_available() -> bool:
    try:
        require_core()
        return True
    except CoreUnavailable:
        if os.environ.get("SPLAS_REQUIRE_CORE") == "1":
            raise
        return False


@unittest.skipUnless(core_available(), "Fortran core is not built")
class CoreTests(unittest.TestCase):
    def short_config(self, closed=False):
        source = (
            ROOT / "configs" / "reference" / "reference_closed_constant_default.json"
            if closed
            else ROOT / "configs" / "reference" / "smoke.json"
        )
        config = load_config(source)
        config["grid"]["layers"] = 7
        config["numerics"].update(
            {
                "duration_days": 1.0,
                "timestep_days": 1.0 / 24.0,
                "output_interval_days": 1.0 / 6.0,
                "seed": 11,
            }
        )
        return config

    def test_deterministic_finite_nonnegative(self):
        config = self.short_config()
        first = run_core(config)
        second = run_core(config)
        for name in ("P", "Z", "quota_P", "nutrient", "detritus", "irradiance"):
            self.assertTrue(np.all(np.isfinite(first[name])))
            self.assertTrue(np.all(first[name] >= 0.0))
            np.testing.assert_array_equal(first[name], second[name])

    def test_closed_mass_balance(self):
        result = run_core(self.short_config(closed=True))
        tolerance = 1e-6 * max(abs(result["mass_balance_error"][0]), 1.0)
        self.assertLessEqual(np.max(np.abs(result["mass_balance_error"])), tolerance)

    def test_timestep_convergence(self):
        coarse = self.short_config()
        fine = copy.deepcopy(coarse)
        coarse["numerics"]["timestep_days"] = 1.0 / 24.0
        fine["numerics"]["timestep_days"] = 1.0 / 48.0
        coarse_result = run_core(coarse)
        fine_result = run_core(fine)
        p_coarse = coarse_result["P"][-1].sum()
        p_fine = fine_result["P"][-1].sum()
        self.assertLess(abs(p_coarse - p_fine) / max(abs(p_fine), 1e-12), 0.05)

    def test_output_frequency_invariance(self):
        first = self.short_config()
        second = copy.deepcopy(first)
        second["numerics"]["output_interval_days"] = 1.0 / 12.0
        a = run_core(first)
        b = run_core(second)
        np.testing.assert_allclose(a["P"][-1], b["P"][-1], rtol=1e-12, atol=1e-14)
