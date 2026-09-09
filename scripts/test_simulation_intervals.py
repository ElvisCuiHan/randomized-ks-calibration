import unittest

from scipy.stats import binomtest

from ks_scdesign_extended_simulations import rate_interval
from refresh_simulation_intervals import reconstruct_count


class IntervalTests(unittest.TestCase):
    def test_wilson_matches_scipy_including_boundaries(self):
        for k, n in ((0, 1000), (1000, 1000), (25, 400), (12, 300), (663, 1000)):
            expected = binomtest(k, n).proportion_ci(method="wilson")
            actual = rate_interval(k / n, n)
            self.assertAlmostEqual(actual[0], expected.low)
            self.assertAlmostEqual(actual[1], expected.high)
            self.assertGreater(actual[1] - actual[0], 0)

    def test_count_recovery(self):
        self.assertEqual(reconstruct_count("0.062 [0.039; 0.086]", 400), 25)
        self.assertEqual(reconstruct_count("0.003 [0.000; 0.010]", 300), 1)
        with self.assertRaises(ValueError):
            reconstruct_count("0.061 [0; 1]", 400)


if __name__ == "__main__":
    unittest.main()
