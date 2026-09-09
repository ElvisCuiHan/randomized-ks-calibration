#!/usr/bin/env python3
"""Extended simulations for randomized KS calibration.

The script produces small reproducible tables for the single-cell simulator
section of the manuscript:

1. Known versus fitted negative-binomial margins across sample sizes, including
   fixed-parameter and refitted bootstrap references and zero-inflation power.
2. A streaming scalability check through one million cells.
3. A copula-only misspecification experiment showing that marginal KS
   diagnostics stay near size while a residual correlation diagnostic reacts.
4. A pointer to the separately refitted conditional experiment implemented in R.
5. A lightweight real count-data application using R's InsectSprays data values
   embedded below, avoiding external data downloads.
"""

from __future__ import annotations

import argparse
import math
import time
from dataclasses import dataclass

import numpy as np

try:
    from scipy.special import ndtr as _normal_cdf
except Exception:  # pragma: no cover - fallback for minimal Python installs.
    _erf = np.vectorize(math.erf)

    def _normal_cdf(x: np.ndarray) -> np.ndarray:
        return 0.5 * (1.0 + _erf(x / math.sqrt(2.0)))


def ks_stat_uniform(u: np.ndarray) -> float:
    u = np.sort(np.asarray(u, dtype=float))
    n = u.size
    upper = np.arange(1, n + 1, dtype=float) / n
    lower = np.arange(0, n, dtype=float) / n
    return float(max(np.max(upper - u), np.max(u - lower)))


def two_sample_ks(x: np.ndarray, y: np.ndarray) -> float:
    x = np.sort(np.asarray(x, dtype=float))
    y = np.sort(np.asarray(y, dtype=float))
    values = np.sort(np.concatenate([x, y]))
    fx = np.searchsorted(x, values, side="right") / x.size
    fy = np.searchsorted(y, values, side="right") / y.size
    return float(np.max(np.abs(fx - fy)))


def two_sample_cutoff(n: int, m: int, alpha: float) -> float:
    return math.sqrt(-0.5 * math.log(alpha / 2.0)) * math.sqrt((n + m) / (n * m))


def kolmogorov_pvalue(d: float, n: int, terms: int = 80) -> float:
    z = math.sqrt(n) * d
    total = 0.0
    for k in range(1, terms + 1):
        total += (-1) ** (k - 1) * math.exp(-2.0 * k * k * z * z)
    return max(0.0, min(1.0, 2.0 * total))


def two_sample_pvalue(d: float, n: int, m: int, terms: int = 80) -> float:
    z = math.sqrt(n * m / (n + m)) * d
    total = 0.0
    for k in range(1, terms + 1):
        total += (-1) ** (k - 1) * math.exp(-2.0 * k * k * z * z)
    return max(0.0, min(1.0, 2.0 * total))


def rate_interval(rate: float, trials: int) -> tuple[float, float]:
    if trials <= 0 or not 0.0 <= rate <= 1.0:
        raise ValueError("A rate in [0, 1] and positive trial count are required")
    z = 1.959963984540054
    denominator = 1.0 + z * z / trials
    center = (rate + z * z / (2.0 * trials)) / denominator
    radius = z * math.sqrt(rate * (1.0 - rate) / trials + z * z / (4.0 * trials**2)) / denominator
    return max(0.0, center - radius), min(1.0, center + radius)


def format_rate(rate: float, trials: int) -> str:
    lo, hi = rate_interval(rate, trials)
    return f"{rate:.3f} [{lo:.3f}; {hi:.3f}]"


def poisson_pit_scalar(count: int, v: float, mean: float) -> float:
    pmf = math.exp(-mean)
    left = 0.0
    for k in range(count):
        left += pmf
        pmf *= mean / (k + 1.0)
    return min(max(left + v * pmf, 1e-12), 1.0 - 1e-12)


def nb_p(size: float, mean: float) -> float:
    return size / (size + mean)


def nb_pmf_table(size: float, mean: float, max_k: int) -> np.ndarray:
    p = nb_p(size, mean)
    pmf = np.empty(max_k + 1, dtype=float)
    pmf[0] = p**size
    for k in range(max_k):
        pmf[k + 1] = pmf[k] * (k + size) / (k + 1.0) * (1.0 - p)
    total = pmf.sum()
    if total > 1.0:
        pmf /= total
    return pmf


def nb_cdf_tables(size: float, mean: float, max_k: int) -> tuple[np.ndarray, np.ndarray]:
    pmf = nb_pmf_table(size, mean, max_k)
    cdf = np.cumsum(pmf)
    return pmf, cdf


def nb_pit(counts: np.ndarray, v: np.ndarray, size: float, mean: float) -> np.ndarray:
    counts = np.asarray(counts, dtype=int)
    max_k = int(max(counts.max(), math.ceil(mean + 12.0 * math.sqrt(mean + mean * mean / size))))
    pmf, cdf = nb_cdf_tables(size, mean, max_k)
    clipped = np.minimum(counts, max_k)
    left = np.where(clipped > 0, cdf[clipped - 1], 0.0)
    values = left + v * pmf[clipped]
    return np.clip(values, 1e-12, 1.0 - 1e-12)


def nb_pit_varying(
    counts: np.ndarray, v: np.ndarray, size: float, means: np.ndarray
) -> np.ndarray:
    return np.array(
        [nb_pit(np.array([x]), np.array([vv]), size, mu)[0] for x, vv, mu in zip(counts, v, means)]
    )


def nb_ppf_from_uniform(u: np.ndarray, size: float, mean: float) -> np.ndarray:
    u = np.asarray(u, dtype=float)
    max_k = max(50, int(math.ceil(mean + 12.0 * math.sqrt(mean + mean * mean / size))))
    while True:
        _, cdf = nb_cdf_tables(size, mean, max_k)
        if cdf[-1] >= max(float(u.max()), 0.999999) or max_k > 5000:
            return np.searchsorted(cdf, u, side="left")
        max_k *= 2


def nb_rvs(rng: np.random.Generator, n: int, size: float, mean: float) -> np.ndarray:
    p = nb_p(size, mean)
    if size > 1e5:
        return rng.poisson(mean, size=n)
    return rng.negative_binomial(size, p, size=n)


def fit_nb_mom(counts: np.ndarray) -> tuple[float, float]:
    counts = np.asarray(counts, dtype=float)
    mean = max(float(counts.mean()), 1e-6)
    var = float(counts.var(ddof=1)) if counts.size > 1 else mean
    if var <= mean * 1.02:
        size = 1000.0
    else:
        size = min(1000.0, max(0.05, mean * mean / (var - mean)))
    return size, mean


def gaussian_copula_nb(
    rng: np.random.Generator,
    n: int,
    rho: float,
    size1: float,
    mean1: float,
    size2: float,
    mean2: float,
) -> tuple[np.ndarray, np.ndarray]:
    z1 = rng.normal(size=n)
    z2 = rho * z1 + math.sqrt(max(0.0, 1.0 - rho * rho)) * rng.normal(size=n)
    u1 = np.clip(_normal_cdf(z1), 1e-12, 1.0 - 1e-12)
    u2 = np.clip(_normal_cdf(z2), 1e-12, 1.0 - 1e-12)
    return nb_ppf_from_uniform(u1, size1, mean1), nb_ppf_from_uniform(u2, size2, mean2)


@dataclass(frozen=True)
class MarginResult:
    model: str
    brownian_cutoff: float
    refit_bootstrap: float | None


@dataclass(frozen=True)
class MarginScaleResult:
    n: int
    null_mc: int
    power_mc: int
    boot: int
    known_margin: float
    naive_ks: float
    fixed_parameter_bootstrap: float
    refitted_bootstrap: float
    misspecification_power: float
    elapsed_seconds: float


@dataclass(frozen=True)
class ScalabilityResult:
    n: int
    genes: int
    observed_seconds: float
    one_refit_seconds: float
    projected_199_serial_minutes: float
    streaming_working_arrays_mb: float


def simulate_fitted_nb_margin(
    rng: np.random.Generator,
    n: int,
    mc: int,
    boot: int,
    alpha: float,
    size: float = 2.5,
    mean: float = 6.0,
) -> list[MarginResult]:
    cutoff = 1.358 / math.sqrt(n)
    known_reject = []
    fitted_reject = []
    bootstrap_reject = []
    for _ in range(mc):
        x = nb_rvs(rng, n, size, mean)
        known_u = nb_pit(x, rng.random(n), size, mean)
        known_reject.append(ks_stat_uniform(known_u) > cutoff)

        fit_size, fit_mean = fit_nb_mom(x)
        fitted_u = nb_pit(x, rng.random(n), fit_size, fit_mean)
        observed = ks_stat_uniform(fitted_u)
        fitted_reject.append(observed > cutoff)

        boot_stats = []
        for _b in range(boot):
            xb = nb_rvs(rng, n, fit_size, fit_mean)
            b_size, b_mean = fit_nb_mom(xb)
            ub = nb_pit(xb, rng.random(n), b_size, b_mean)
            boot_stats.append(ks_stat_uniform(ub))
        boot_cutoff = float(np.quantile(boot_stats, 1.0 - alpha, method="higher"))
        bootstrap_reject.append(observed > boot_cutoff)

    return [
        MarginResult("Known NB margin", float(np.mean(known_reject)), None),
        MarginResult("Fitted NB margin", float(np.mean(fitted_reject)), float(np.mean(bootstrap_reject))),
    ]


def zero_inflated_nb_rvs(
    rng: np.random.Generator,
    n: int,
    size: float,
    mean: float,
    zero_probability: float = 0.15,
) -> np.ndarray:
    counts = nb_rvs(rng, n, size, mean)
    counts[rng.random(n) < zero_probability] = 0
    return counts


def bootstrap_margin_cutoffs(
    rng: np.random.Generator,
    n: int,
    boot: int,
    alpha: float,
    fit_size: float,
    fit_mean: float,
) -> tuple[float, float]:
    fixed_stats = np.empty(boot, dtype=float)
    refit_stats = np.empty(boot, dtype=float)
    for b in range(boot):
        xb = nb_rvs(rng, n, fit_size, fit_mean)
        fixed_u = nb_pit(xb, rng.random(n), fit_size, fit_mean)
        fixed_stats[b] = ks_stat_uniform(fixed_u)

        b_size, b_mean = fit_nb_mom(xb)
        refit_u = nb_pit(xb, rng.random(n), b_size, b_mean)
        refit_stats[b] = ks_stat_uniform(refit_u)

    return (
        float(np.quantile(fixed_stats, 1.0 - alpha, method="higher")),
        float(np.quantile(refit_stats, 1.0 - alpha, method="higher")),
    )


def simulate_margin_scale(
    rng: np.random.Generator,
    n: int,
    null_mc: int,
    power_mc: int,
    boot: int,
    alpha: float,
    size: float = 2.5,
    mean: float = 6.0,
) -> MarginScaleResult:
    start = time.perf_counter()
    cutoff = 1.358 / math.sqrt(n)
    known_reject = np.empty(null_mc, dtype=bool)
    naive_reject = np.empty(null_mc, dtype=bool)
    fixed_boot_reject = np.empty(null_mc, dtype=bool)
    refit_reject = np.empty(null_mc, dtype=bool)

    for trial in range(null_mc):
        x = nb_rvs(rng, n, size, mean)
        known_u = nb_pit(x, rng.random(n), size, mean)
        known_reject[trial] = ks_stat_uniform(known_u) > cutoff

        fit_size, fit_mean = fit_nb_mom(x)
        observed = ks_stat_uniform(nb_pit(x, rng.random(n), fit_size, fit_mean))
        naive_reject[trial] = observed > cutoff
        fixed_cutoff, refit_cutoff = bootstrap_margin_cutoffs(
            rng, n, boot, alpha, fit_size, fit_mean
        )
        fixed_boot_reject[trial] = observed > fixed_cutoff
        refit_reject[trial] = observed > refit_cutoff

    power_reject = np.empty(power_mc, dtype=bool)
    for trial in range(power_mc):
        x = zero_inflated_nb_rvs(rng, n, size, mean)
        fit_size, fit_mean = fit_nb_mom(x)
        observed = ks_stat_uniform(nb_pit(x, rng.random(n), fit_size, fit_mean))
        _, refit_cutoff = bootstrap_margin_cutoffs(
            rng, n, boot, alpha, fit_size, fit_mean
        )
        power_reject[trial] = observed > refit_cutoff

    return MarginScaleResult(
        n=n,
        null_mc=null_mc,
        power_mc=power_mc,
        boot=boot,
        known_margin=float(np.mean(known_reject)),
        naive_ks=float(np.mean(naive_reject)),
        fixed_parameter_bootstrap=float(np.mean(fixed_boot_reject)),
        refitted_bootstrap=float(np.mean(refit_reject)),
        misspecification_power=float(np.mean(power_reject)),
        elapsed_seconds=time.perf_counter() - start,
    )


def benchmark_streaming_scalability(
    rng: np.random.Generator,
    sizes: list[int],
    genes: int,
    size: float = 2.5,
    mean: float = 6.0,
) -> list[ScalabilityResult]:
    results: list[ScalabilityResult] = []
    for n in sizes:
        fitted_parameters: list[tuple[float, float]] = []
        start = time.perf_counter()
        for _gene in range(genes):
            x = nb_rvs(rng, n, size, mean)
            fit_size, fit_mean = fit_nb_mom(x)
            fitted_parameters.append((fit_size, fit_mean))
            ks_stat_uniform(nb_pit(x, rng.random(n), fit_size, fit_mean))
        observed_seconds = time.perf_counter() - start

        start = time.perf_counter()
        for fit_size, fit_mean in fitted_parameters:
            xb = nb_rvs(rng, n, fit_size, fit_mean)
            b_size, b_mean = fit_nb_mom(xb)
            ks_stat_uniform(nb_pit(xb, rng.random(n), b_size, b_mean))
        one_refit_seconds = time.perf_counter() - start

        results.append(
            ScalabilityResult(
                n=n,
                genes=genes,
                observed_seconds=observed_seconds,
                one_refit_seconds=one_refit_seconds,
                projected_199_serial_minutes=(observed_seconds + 199 * one_refit_seconds) / 60.0,
                streaming_working_arrays_mb=24.0 * n / (1024.0 * 1024.0),
            )
        )
    return results


@dataclass(frozen=True)
class CopulaResult:
    scenario: str
    gene1_ks: float
    gene2_ks: float
    residual_corr: float


def simulate_copula_blind_spot(
    rng: np.random.Generator,
    n: int,
    mc: int,
    alpha: float,
    rho_real: float = 0.65,
) -> list[CopulaResult]:
    scenarios = [
        ("same copula (rho=0.65)", rho_real),
        ("moderate wrong copula (rho=0.50)", 0.50),
        ("rho=0 independence copula", 0.0),
    ]
    out: list[CopulaResult] = []
    cutoff = two_sample_cutoff(n, n, alpha)
    se = math.sqrt(2.0 / (n - 3.0))
    for label, rho_sim in scenarios:
        gene1_reject = []
        gene2_reject = []
        corr_reject = []
        for _ in range(mc):
            x1, x2 = gaussian_copula_nb(rng, n, rho_real, 2.0, 5.0, 3.0, 8.0)
            y1, y2 = gaussian_copula_nb(rng, n, rho_sim, 2.0, 5.0, 3.0, 8.0)
            ux1 = nb_pit(x1, rng.random(n), 2.0, 5.0)
            ux2 = nb_pit(x2, rng.random(n), 3.0, 8.0)
            uy1 = nb_pit(y1, rng.random(n), 2.0, 5.0)
            uy2 = nb_pit(y2, rng.random(n), 3.0, 8.0)

            gene1_reject.append(two_sample_ks(ux1, uy1) > cutoff)
            gene2_reject.append(two_sample_ks(ux2, uy2) > cutoff)

            rx = np.corrcoef(ux1, ux2)[0, 1]
            ry = np.corrcoef(uy1, uy2)[0, 1]
            rx = float(np.clip(rx, -0.999, 0.999))
            ry = float(np.clip(ry, -0.999, 0.999))
            z = abs(math.atanh(rx) - math.atanh(ry)) / se
            corr_reject.append(z > 1.96)
        out.append(
            CopulaResult(
                label,
                float(np.mean(gene1_reject)),
                float(np.mean(gene2_reject)),
                float(np.mean(corr_reject)),
            )
        )
    return out


def insect_sprays_application(seed: int) -> tuple[int, float, float, float]:
    counts = np.array(
        [
            10, 7, 20, 14, 14, 12, 10, 23, 17, 20, 14, 13,
            11, 17, 21, 11, 16, 14, 17, 17, 19, 21, 7, 13,
            0, 1, 7, 2, 3, 1, 2, 1, 3, 0, 1, 4,
            3, 5, 12, 6, 4, 3, 5, 5, 5, 5, 2, 4,
            3, 5, 3, 5, 3, 6, 1, 1, 3, 2, 6, 4,
            11, 9, 15, 22, 15, 16, 13, 10, 26, 26, 24, 13,
        ],
        dtype=int,
    )
    sprays = np.repeat(np.arange(6), 12)
    rng = np.random.default_rng(seed)
    means = np.array([counts[sprays == g].mean() for g in range(6)])
    fitted = means[sprays]
    u = np.array(
        [
            poisson_pit_scalar(int(x), float(v), float(mu))
            for x, v, mu in zip(counts, rng.random(counts.size), fitted)
        ]
    )
    d = ks_stat_uniform(u)
    pvalue = kolmogorov_pvalue(d, counts.size)
    max_group = 0.0
    for i in range(6):
        for j in range(i + 1, 6):
            max_group = max(max_group, two_sample_ks(u[sprays == i], u[sprays == j]))
    return counts.size, d, pvalue, max_group


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=120)
    parser.add_argument("--mc", type=int, default=1000)
    parser.add_argument("--boot", type=int, default=499)
    parser.add_argument("--scale-sizes", default="250,1000,10000")
    parser.add_argument("--scale-null-mc", default="400,250,300")
    parser.add_argument("--scale-power-mc", default="200,125,100")
    parser.add_argument("--scale-boot", default="199,149,99")
    parser.add_argument("--scalability-sizes", default="1000,10000,100000,1000000")
    parser.add_argument("--scalability-genes", type=int, default=10)
    parser.add_argument(
        "--skip-revision-scale",
        action="store_true",
        help="Skip the expanded size, misspecification, and scalability experiments.",
    )
    parser.add_argument("--alpha", type=float, default=0.05)
    parser.add_argument("--seed", type=int, default=20260609)
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)

    if not args.skip_revision_scale:
        scale_sizes = [int(value) for value in args.scale_sizes.split(",")]
        scale_null_mc = [int(value) for value in args.scale_null_mc.split(",")]
        scale_power_mc = [int(value) for value in args.scale_power_mc.split(",")]
        scale_boot = [int(value) for value in args.scale_boot.split(",")]
        lengths = {len(scale_sizes), len(scale_null_mc), len(scale_power_mc), len(scale_boot)}
        if lengths != {len(scale_sizes)}:
            raise SystemExit("All --scale-* comma-separated lists must have equal length.")

        print("fitted_nb_scale")
        print(
            "n,null_mc,power_mc,boot,known_margin_with_95mc_ci,"
            "naive_ks_with_95mc_ci,fixed_parameter_bootstrap_with_95mc_ci,"
            "refitted_bootstrap_with_95mc_ci,zero_inflated_power_with_95mc_ci,"
            "elapsed_seconds"
        )
        for n_scale, null_mc, power_mc, boot_scale in zip(
            scale_sizes, scale_null_mc, scale_power_mc, scale_boot
        ):
            result = simulate_margin_scale(
                rng,
                n=n_scale,
                null_mc=null_mc,
                power_mc=power_mc,
                boot=boot_scale,
                alpha=args.alpha,
            )
            print(
                f"{result.n},{result.null_mc},{result.power_mc},{result.boot},"
                f"{format_rate(result.known_margin, result.null_mc)},"
                f"{format_rate(result.naive_ks, result.null_mc)},"
                f"{format_rate(result.fixed_parameter_bootstrap, result.null_mc)},"
                f"{format_rate(result.refitted_bootstrap, result.null_mc)},"
                f"{format_rate(result.misspecification_power, result.power_mc)},"
                f"{result.elapsed_seconds:.3f}"
            )

        print("\nlarge_scale_feasibility")
        print(
            "n,genes,observed_fit_residual_seconds,one_refit_seconds,"
            "projected_199_serial_minutes,streaming_working_arrays_mb"
        )
        scalability_sizes = [int(value) for value in args.scalability_sizes.split(",")]
        for result in benchmark_streaming_scalability(
            rng, scalability_sizes, args.scalability_genes
        ):
            print(
                f"{result.n},{result.genes},{result.observed_seconds:.4f},"
                f"{result.one_refit_seconds:.4f},"
                f"{result.projected_199_serial_minutes:.3f},"
                f"{result.streaming_working_arrays_mb:.2f}"
            )

    print(f"\nsettings,n={args.n},mc={args.mc},boot={args.boot},alpha={args.alpha},seed={args.seed}")
    print("\nfitted_nb_margin")
    print("model,brownian_cutoff_with_95mc_ci,refit_bootstrap_with_95mc_ci")
    for result in simulate_fitted_nb_margin(rng, args.n, args.mc, args.boot, args.alpha):
        bootstrap = "--" if result.refit_bootstrap is None else format_rate(result.refit_bootstrap, args.mc)
        print(f"{result.model},{format_rate(result.brownian_cutoff, args.mc)},{bootstrap}")

    print("\ncopula_blind_spot")
    print("scenario,gene1_marginal_ks_with_95mc_ci,gene2_marginal_ks_with_95mc_ci,residual_correlation_with_95mc_ci")
    for result in simulate_copula_blind_spot(rng, 2 * args.n, args.mc, args.alpha):
        print(
            f"{result.scenario},{format_rate(result.gene1_ks, args.mc)},"
            f"{format_rate(result.gene2_ks, args.mc)},{format_rate(result.residual_corr, args.mc)}"
        )

    print("\nconditional_calibration_source")
    print("Run scripts/ks_conditional_calibration.R; see results/conditional_calibration/summary.csv")

    print("\nreal_count_application")
    print("dataset,n,ks_stat,uniform_pvalue,max_pairwise_spray_ks")
    n, d, pvalue, max_group = insect_sprays_application(args.seed + 17)
    print(f"InsectSprays,{n},{d:.3f},{pvalue:.3f},{max_group:.3f}")


if __name__ == "__main__":
    main()
