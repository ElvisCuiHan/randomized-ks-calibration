"""Focused numerical and reproducibility checks for the YL revision."""

import tempfile
import unittest
from pathlib import Path

import numpy as np

from imvigor210_residual_diagnostics import (
    analyze_gene, fit_conditional_nb, holm_adjust, monte_carlo_interval,
)


class BootstrapTests(unittest.TestCase):
    def setUp(self):
        self.groups = np.repeat(np.arange(2), 100)
        self.exposure = np.linspace(0.5, 2, 200)
        self.counts = np.random.default_rng(2).negative_binomial(
            2, 2 / (2 + 4 * self.exposure)
        )

    def test_shared_moment_formula(self):
        fit = fit_conditional_nb(self.counts, self.exposure, self.groups, 2)
        rates = np.array([
            self.counts[self.groups == j].sum() / self.exposure[self.groups == j].sum()
            for j in range(2)
        ])
        means = self.exposure * rates[self.groups]
        alpha = max(np.sum((self.counts - means) ** 2 - self.counts) / np.sum(means**2), .001)
        np.testing.assert_allclose(fit.means, means)
        self.assertAlmostEqual(fit.size, np.clip(1 / alpha, .05, 1000))

    def test_patient_dispersion_equals_separate_fits(self):
        fit = fit_conditional_nb(self.counts, self.exposure, self.groups, 2, "patient")
        for j in range(2):
            mask = self.groups == j
            separate = fit_conditional_nb(self.counts[mask], self.exposure[mask],
                                          np.zeros(mask.sum(), dtype=int), 1)
            np.testing.assert_allclose(fit.size[mask], separate.size)

    def test_observed_pit_and_bootstrap_prefix(self):
        with tempfile.TemporaryDirectory() as tmp:
            for model in ("shared", "patient"):
                first = analyze_gene("test", self.counts, self.exposure, self.groups,
                                     2, 9, 42, model, tmp + "/a")
                second = analyze_gene("test", self.counts, self.exposure, self.groups,
                                      2, 19, 42, model, tmp + "/b")
                a = np.load(Path(tmp) / "a/replicates" / f"{model}_test.npz")
                b = np.load(Path(tmp) / "b/replicates" / f"{model}_test.npz")
                for key in ("fixed_statistics", "refitted_statistics"):
                    np.testing.assert_array_equal(a[key], b[key][:9])
                np.testing.assert_array_equal(a["observed_randomized_residuals"],
                                              b["observed_randomized_residuals"])
                self.assertEqual(first["observed_stratified_ks"], second["observed_stratified_ks"])
                self.assertEqual(first["refitted_p"], (1 + first["refitted_exceedances"]) / 10)

    def test_holm_and_mc_boundaries(self):
        np.testing.assert_allclose(holm_adjust(np.array([.001, .02, .2])), [.003, .04, .2])
        self.assertEqual(monte_carlo_interval(0, 9999)[0], 0)
        self.assertLess(monte_carlo_interval(0, 9999)[1], .0004)
        self.assertEqual(monte_carlo_interval(9999, 9999)[1], 1)


if __name__ == "__main__":
    unittest.main()
