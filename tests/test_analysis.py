import unittest

import numpy as np

from splas.analysis import (
    carbon_from_phyto,
    critical_depth,
    detect_blooms,
    euphotic_depth,
    export_efficiencies,
    shannon_metrics,
)


class AnalysisTests(unittest.TestCase):
    def test_shannon_uniform_and_monodominant(self):
        uniform = np.ones((1, 4))
        entropy, evenness, count, valid = shannon_metrics(uniform)
        self.assertAlmostEqual(entropy[0], np.log(4.0))
        self.assertAlmostEqual(evenness[0], 1.0)
        self.assertEqual(count[0], 4)
        self.assertTrue(valid[0])
        mono = np.array([[1.0, 0.0, 0.0, 0.0]])
        entropy, evenness, count, valid = shannon_metrics(mono)
        self.assertAlmostEqual(entropy[0], 0.0)
        self.assertAlmostEqual(evenness[0], 0.0)
        self.assertEqual(count[0], 1)

    def test_low_biomass_diversity_is_masked(self):
        entropy, evenness, count, valid = shannon_metrics(
            np.full((1, 4), 1e-20), community_threshold=1e-10
        )
        self.assertFalse(valid[0])
        self.assertTrue(np.isnan(entropy[0]))
        self.assertTrue(np.isnan(evenness[0]))

    def test_euphotic_depth(self):
        depth = np.array([2.5, 7.5, 12.5, 17.5])
        light = np.array([[100.0, 20.0, 1.1, 0.9]])
        result = euphotic_depth(light, np.array([100.0]), depth)
        self.assertAlmostEqual(result[0], 12.5)

    def test_critical_depth_synthetic_crossing(self):
        time = np.array([0.0, 1.0])
        depth = np.array([5.0, 10.0, 15.0])
        initial = np.array([1.0, 1.0, 10.0])
        final = np.array([2.0, 1.0, 0.1])
        crossing, growth = critical_depth(np.vstack([initial, final]), time, depth)
        self.assertTrue(np.isfinite(crossing[1]))
        self.assertGreaterEqual(crossing[1], depth[0])
        self.assertLessEqual(crossing[1], depth[-1])

    def test_carbon_units_and_ratio(self):
        p = np.array([[[1.0]]])
        quota = np.array([[[0.5]]])
        carbon = np.array([12.0])
        cells, carbon_mmol, ratio = carbon_from_phyto(p, quota, carbon)
        self.assertAlmostEqual(cells.item(), 2.0)
        self.assertAlmostEqual(carbon_mmol.item(), 2e-9)
        self.assertAlmostEqual(ratio.item(), 2e-9)

    def test_export_metrics_simple_case(self):
        time = np.array([0.0, 1.0, 2.0])
        depth = np.array([0.0, 1.0])
        export = np.array([2.0, 2.0, 2.0])
        biomass = np.ones((3, 2))
        light = np.full(3, 4.0)
        metrics = export_efficiencies(export, biomass, light, time, depth, transient_days=0.0)
        self.assertAlmostEqual(metrics["export_mean_mg_C_m2_d"], 2.0)
        self.assertAlmostEqual(metrics["export_integral_mg_C_m2"], 4.0)
        self.assertAlmostEqual(metrics["export_per_biomass"], 2.0)
        self.assertAlmostEqual(metrics["export_per_light"], 0.5)

    def test_bloom_onset_definition(self):
        time = np.arange(6.0)
        depth = np.array([0.0, 1.0])
        integrated = np.array([2.0, 1.0, 2.0, 4.0, 2.0, 1.0])
        biomass = np.column_stack([integrated / 2.0, integrated / 2.0])
        events = detect_blooms(biomass, time, depth, minimum_biomass=0.1)
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0]["onset_day"], 2.0)
        self.assertEqual(events[0]["peak_day"], 3.0)
