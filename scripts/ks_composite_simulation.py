#!/usr/bin/env python3
"""Monte Carlo calibration checks for the composite-null KS discussion.

The script reproduces Table 1 in the manuscript.  It uses synthetic data only:
known standard normal samples, normal samples with parameters re-estimated, and
exponential samples with the rate re-estimated.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass

import numpy as np

try:
    from scipy.special import ndtr as _normal_cdf
except Exception:  # pragma: no cover - fallback for minimal Python installs.
    _erf = np.vectorize(math.erf)

    def _normal_cdf(x: np.ndarray) -> np.ndarray:
        return 0.5 * (1.0 + _erf(x / math.sqrt(2.0)))


@dataclass(frozen=True)
class Result:
    model: str
    classical: float
    bootstrap: float | None


def ks_stat_from_cdf(cdf_values: np.ndarray) -> np.ndarray:
    """Two-sided one-sample KS statistic for sorted row-wise CDF values."""
    if cdf_values.ndim == 1:
        cdf_values = cdf_values[None, :]
    n = cdf_values.shape[1]
    upper = np.arange(1, n + 1, dtype=float) / n
    lower = np.arange(0, n, dtype=float) / n
    d_plus = np.max(upper - cdf_values, axis=1)
    d_minus = np.max(cdf_values - lower, axis=1)
    return np.maximum(d_plus, d_minus)


def ks_known_normal(samples: np.ndarray) -> np.ndarray:
    xs = np.sort(samples, axis=1)
    return ks_stat_from_cdf(_normal_cdf(xs))


def ks_fitted_normal(samples: np.ndarray) -> np.ndarray:
    xs = np.sort(samples, axis=1)
    mu = samples.mean(axis=1, keepdims=True)
    sigma = samples.std(axis=1, ddof=0, keepdims=True)
    return ks_stat_from_cdf(_normal_cdf((xs - mu) / sigma))


def ks_fitted_exponential(samples: np.ndarray) -> np.ndarray:
    xs = np.sort(samples, axis=1)
    rate = 1.0 / samples.mean(axis=1, keepdims=True)
    cdf = 1.0 - np.exp(-rate * xs)
    return ks_stat_from_cdf(cdf)


def bootstrap_rejection_normal(
    rng: np.random.Generator,
    sample: np.ndarray,
    observed: float,
    boot: int,
    alpha: float,
) -> bool:
    mu = sample.mean()
    sigma = sample.std(ddof=0)
    boot_samples = rng.normal(mu, sigma, size=(boot, sample.size))
    cutoff = np.quantile(ks_fitted_normal(boot_samples), 1.0 - alpha, method="higher")
    return bool(observed > cutoff)


def bootstrap_rejection_exponential(
    rng: np.random.Generator,
    sample: np.ndarray,
    observed: float,
    boot: int,
    alpha: float,
) -> bool:
    scale = sample.mean()
    boot_samples = rng.exponential(scale, size=(boot, sample.size))
    cutoff = np.quantile(ks_fitted_exponential(boot_samples), 1.0 - alpha, method="higher")
    return bool(observed > cutoff)


def run_simulation(n: int, mc: int, boot: int, alpha: float, seed: int) -> list[Result]:
    rng = np.random.default_rng(seed)
    classical_cutoff = 1.358 / math.sqrt(n)

    known = rng.normal(size=(mc, n))
    known_rate = float(np.mean(ks_known_normal(known) > classical_cutoff))

    normal = rng.normal(size=(mc, n))
    normal_stats = ks_fitted_normal(normal)
    normal_classical = float(np.mean(normal_stats > classical_cutoff))
    normal_boot = np.mean(
        [
            bootstrap_rejection_normal(rng, normal[i], normal_stats[i], boot, alpha)
            for i in range(mc)
        ]
    )

    exponential = rng.exponential(1.0, size=(mc, n))
    exponential_stats = ks_fitted_exponential(exponential)
    exponential_classical = float(np.mean(exponential_stats > classical_cutoff))
    exponential_boot = np.mean(
        [
            bootstrap_rejection_exponential(rng, exponential[i], exponential_stats[i], boot, alpha)
            for i in range(mc)
        ]
    )

    return [
        Result("Known N(0,1)", known_rate, None),
        Result("Fitted normal (mu, sigma^2)", normal_classical, float(normal_boot)),
        Result("Fitted exponential rate", exponential_classical, float(exponential_boot)),
    ]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=100)
    parser.add_argument("--mc", type=int, default=2000)
    parser.add_argument("--boot", type=int, default=399)
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=20260520)
    args = parser.parse_args()

    results = run_simulation(args.n, args.mc, args.boot, args.alpha, args.seed)
    print(f"n={args.n}, mc={args.mc}, boot={args.boot}, alpha={args.alpha}, seed={args.seed}")
    print("model,classical,bootstrap")
    for result in results:
        bootstrap = "--" if result.bootstrap is None else f"{result.bootstrap:.3f}"
        print(f"{result.model},{result.classical:.3f},{bootstrap}")


if __name__ == "__main__":
    main()
