import unittest

import numpy as np

from splas.forcing import attenuated_irradiance, light_limitation, surface_irradiance
from splas.traits import (
    droop,
    log_size_grid,
    monod,
    phyto_traits,
    preference_matrix,
    zoo_traits,
)


class TraitTests(unittest.TestCase):
    def test_log_spaced_size_ranges(self):
        phyto = log_size_grid(1.0, 200.0)
        zoo = log_size_grid(2.0, 1000.0)
        self.assertEqual(len(phyto), 70)
        self.assertEqual(len(zoo), 70)
        self.assertAlmostEqual(phyto[0], 1.0)
        self.assertAlmostEqual(phyto[-1], 200.0)
        self.assertAlmostEqual(zoo[0], 2.0)
        self.assertAlmostEqual(zoo[-1], 1000.0)
        np.testing.assert_allclose(np.diff(np.log10(phyto)), np.diff(np.log10(phyto))[0])

    def test_allometries_and_preference(self):
        phyto = log_size_grid(1.0, 200.0)
        zoo = log_size_grid(2.0, 1000.0)
        ptraits = phyto_traits(phyto)
        ztraits = zoo_traits(zoo)
        preference = preference_matrix(phyto, ztraits["optimal_prey_um"], 0.25)
        self.assertTrue(np.all(np.diff(ptraits["cell_volume_um3"]) > 0))
        self.assertTrue(np.all(np.diff(ptraits["mu_max_d"]) < 0))
        for j in (0, 20, 40, 69):
            preferred = phyto[np.argmax(preference[:, j])]
            nearest = phyto[np.argmin(abs(np.log10(phyto) - np.log10(ztraits["optimal_prey_um"][j])))]
            self.assertEqual(preferred, nearest)

    def test_monod_droop_limits(self):
        half = np.array([1.0, 2.0])
        np.testing.assert_allclose(monod(0.0, half), 0.0)
        np.testing.assert_allclose(monod(1e12, half), 1.0, rtol=1e-10)
        mu = np.array([1.0, 2.0])
        qmin = np.array([0.1, 0.2])
        np.testing.assert_allclose(droop(mu, qmin, qmin), 0.0)
        np.testing.assert_allclose(droop(mu, qmin, qmin * 2), mu / 2)


class ForcingTests(unittest.TestCase):
    def test_constant_and_seasonal_forcing(self):
        time = np.array([0.0, 365.0 / 4.0, 365.0 / 2.0])
        constant = surface_irradiance(time, "constant", reference=100.0)
        seasonal = surface_irradiance(time, "seasonal", reference=100.0)
        np.testing.assert_allclose(constant, 100.0)
        self.assertAlmostEqual(seasonal[0], 50.0)
        self.assertAlmostEqual(seasonal[1], 100.0)
        self.assertAlmostEqual(seasonal[2], 50.0)

    def test_daily_linear_and_normal_are_deterministic(self):
        time = np.linspace(0.0, 2.0, 9)
        daily = surface_irradiance(time, "daily", reference=100.0)
        self.assertAlmostEqual(daily[0], 0.0)
        self.assertAlmostEqual(daily[2], 100.0)
        linear = surface_irradiance(time, "linear", reference=100.0, total_delta=20.0, duration_days=2.0)
        self.assertAlmostEqual(linear[-1], 120.0)
        first = surface_irradiance(time, "normal", reference=100.0, sd=5.0, seed=4)
        second = surface_irradiance(time, "normal", reference=100.0, sd=5.0, seed=4)
        np.testing.assert_array_equal(first, second)

    def test_attenuation_and_self_shading(self):
        depth = np.array([2.5, 7.5, 12.5])
        low = attenuated_irradiance(100.0, depth, np.array([0.1, 0.1, 0.1]), 5.0)
        high = attenuated_irradiance(100.0, depth, np.array([1.0, 1.0, 1.0]), 5.0)
        self.assertTrue(np.all(np.diff(low) < 0))
        self.assertTrue(np.all(high < low))
        limitation = light_limitation(np.array([0.0, 1000.0]), np.array([1.0]))
        self.assertAlmostEqual(limitation[0, 0], 0.0)
        self.assertGreater(limitation[1, 0], 0.99)
