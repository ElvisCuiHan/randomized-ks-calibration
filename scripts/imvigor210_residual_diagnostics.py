#!/usr/bin/env python3
"""Residual calibration for public IMvigor210 baseline PBMC counts.

The analysis uses the per-patient integer count matrices deposited under GEO
GSE145281. It diagnoses a deliberately transparent working model: a
patient-specific negative-binomial mean with observed library size as a fixed
offset and shared or patient-specific gene dispersions. Clinical response labels are used only
to describe the ten samples; they are not tested or interpreted causally.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import time
import urllib.request
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch
from scipy.stats import beta, nbinom


GEO_SERIES = "GSE145281"
TRIAL_ID = "NCT02108652"
PANEL = (
    "CXCL8",
    "B2M",
    "CD74",
    "HLA-A",
    "HLA-C",
    "HLA-DQA1",
    "HLA-DQB1",
    "HLA-DRB1",
    "TAP1",
    "TAP2",
)
MITOCHONDRIAL_GENES = {
    "ATP6",
    "ATP8",
    "COX1",
    "COX2",
    "COX3",
    "CYTB",
    "ND1",
    "ND2",
    "ND3",
    "ND4",
    "ND5",
    "ND6",
    "RNR1",
    "RNR2",
}
SAMPLES = (
    ("GSM4314121", "R1", "Responder"),
    ("GSM4314122", "R2", "Responder"),
    ("GSM4314123", "R3", "Responder"),
    ("GSM4314124", "R4", "Responder"),
    ("GSM4314125", "R5", "Responder"),
    ("GSM4314126", "NR1", "Nonresponder"),
    ("GSM4314127", "NR2", "Nonresponder"),
    ("GSM4314128", "NR3", "Nonresponder"),
    ("GSM4314129", "NR4", "Nonresponder"),
    ("GSM4314130", "NR5", "Nonresponder"),
)


@dataclass(frozen=True)
class SampleData:
    sample: str
    response: str
    library_size: np.ndarray
    detected_genes: np.ndarray
    mitochondrial_fraction: np.ndarray
    panel_counts: dict[str, np.ndarray]


@dataclass(frozen=True)
class ConditionalNBFit:
    group_rates: np.ndarray
    size: float | np.ndarray
    means: np.ndarray


def sample_path(data_dir: Path, accession: str, sample: str) -> Path:
    return data_dir / f"{accession}_{sample}_raw.txt.gz"


def sample_url(accession: str, sample: str) -> str:
    return (
        "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4314nnn/"
        f"{accession}/suppl/{accession}_{sample}_raw.txt.gz"
    )


def sha256sum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download_missing(data_dir: Path) -> None:
    data_dir.mkdir(parents=True, exist_ok=True)
    for accession, sample, _response in SAMPLES:
        destination = sample_path(data_dir, accession, sample)
        if destination.exists() and destination.stat().st_size > 0:
            continue
        print(f"Downloading {sample} from GEO", flush=True)
        urllib.request.urlretrieve(sample_url(accession, sample), destination)
        with gzip.open(destination, "rb") as handle:
            handle.read(1)


def read_sample(path: Path, sample: str, response: str) -> SampleData:
    library_size: np.ndarray | None = None
    detected_genes: np.ndarray | None = None
    mitochondrial_counts: np.ndarray | None = None
    panel_counts: dict[str, np.ndarray] = {}

    for chunk in pd.read_csv(
        path,
        sep="\t",
        index_col=0,
        compression="gzip",
        chunksize=512,
        dtype={0: str},
    ):
        values = chunk.to_numpy(dtype=np.int64, copy=False)
        if library_size is None:
            cells = values.shape[1]
            library_size = np.zeros(cells, dtype=np.int64)
            detected_genes = np.zeros(cells, dtype=np.int32)
            mitochondrial_counts = np.zeros(cells, dtype=np.int64)

        library_size += values.sum(axis=0, dtype=np.int64)
        detected_genes += (values > 0).sum(axis=0, dtype=np.int32)
        gene_names = chunk.index.astype(str)
        mitochondrial = gene_names.str.startswith("MT-") | gene_names.isin(
            MITOCHONDRIAL_GENES
        )
        if np.any(mitochondrial):
            mitochondrial_counts += values[mitochondrial].sum(axis=0, dtype=np.int64)

        for gene in PANEL:
            if gene in chunk.index:
                panel_counts[gene] = chunk.loc[gene].to_numpy(dtype=np.int64)

    if library_size is None or detected_genes is None or mitochondrial_counts is None:
        raise ValueError(f"No count data found in {path}")
    missing = sorted(set(PANEL).difference(panel_counts))
    if missing:
        raise ValueError(f"Missing prespecified genes in {path}: {', '.join(missing)}")

    mitochondrial_fraction = mitochondrial_counts / np.maximum(library_size, 1)
    return SampleData(
        sample=sample,
        response=response,
        library_size=library_size,
        detected_genes=detected_genes,
        mitochondrial_fraction=mitochondrial_fraction,
        panel_counts=panel_counts,
    )


def qc_mask(sample: SampleData) -> np.ndarray:
    return (
        (sample.library_size >= 500)
        & (sample.detected_genes >= 200)
        & (sample.detected_genes <= 6000)
        & (sample.mitochondrial_fraction <= 0.20)
    )


def fit_conditional_nb(
    counts: np.ndarray,
    exposure: np.ndarray,
    groups: np.ndarray,
    n_groups: int,
    dispersion: str = "shared",
) -> ConditionalNBFit:
    if dispersion not in ("shared", "patient"):
        raise ValueError(f"Unknown dispersion model: {dispersion}")
    group_rates = np.bincount(groups, weights=counts, minlength=n_groups) / np.maximum(
        np.bincount(groups, weights=exposure, minlength=n_groups), 1e-12
    )
    group_rates = np.maximum(group_rates, 1e-8)
    means = exposure * group_rates[groups]

    residual_variance = (counts - means) ** 2 - counts
    if dispersion == "shared":
        alpha = max(float(np.sum(residual_variance) / max(np.sum(means**2), 1e-12)), 0.001)
        size = float(np.clip(1.0 / alpha, 0.05, 1000.0))
    else:
        numerator = np.bincount(groups, weights=residual_variance, minlength=n_groups)
        denominator = np.bincount(groups, weights=means**2, minlength=n_groups)
        alpha = np.maximum(numerator / np.maximum(denominator, 1e-12), 0.001)
        size = np.clip(1.0 / alpha, 0.05, 1000.0)[groups]
    return ConditionalNBFit(group_rates=group_rates, size=size, means=means)


def randomized_nb_pit(
    counts: np.ndarray,
    means: np.ndarray,
    size: float | np.ndarray,
    randomizers: np.ndarray,
) -> np.ndarray:
    probability = size / (size + means)
    left = nbinom.cdf(counts - 1, size, probability)
    mass = nbinom.pmf(counts, size, probability)
    return np.clip(left + randomizers * mass, 1e-12, 1.0 - 1e-12)


def ks_stat_uniform(values: np.ndarray) -> float:
    ordered = np.sort(np.asarray(values, dtype=float))
    n = ordered.size
    upper = np.arange(1, n + 1, dtype=float) / n
    lower = np.arange(0, n, dtype=float) / n
    return float(max(np.max(upper - ordered), np.max(ordered - lower)))


def stratified_ks_stat(values: np.ndarray, groups: np.ndarray, n_groups: int) -> float:
    return max(
        np.sqrt(np.sum(groups == group)) * ks_stat_uniform(values[groups == group])
        for group in range(n_groups)
    )


def simulate_nb(
    rng: np.random.Generator,
    means: np.ndarray,
    size: float | np.ndarray,
) -> np.ndarray:
    probability = size / (size + means)
    return rng.negative_binomial(size, probability).astype(np.int64)


def holm_adjust(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values)
    adjusted = np.empty_like(p_values, dtype=float)
    running = 0.0
    total = len(p_values)
    for rank, index in enumerate(order):
        running = max(running, (total - rank) * p_values[index])
        adjusted[index] = min(running, 1.0)
    return adjusted


def analyze_gene(
    gene: str,
    counts: np.ndarray,
    exposure: np.ndarray,
    groups: np.ndarray,
    n_groups: int,
    boot: int,
    seed: int,
    dispersion: str = "shared",
    output_dir: str | None = None,
) -> dict[str, float | int | str]:
    if boot < 1:
        raise ValueError("At least one bootstrap replicate is required")
    # Observed randomizers are independent of B, worker order, and candidate model.
    observed_rng = np.random.default_rng(seed)
    rng = np.random.default_rng(np.random.SeedSequence([seed, 20260907]))
    observed_fit = fit_conditional_nb(counts, exposure, groups, n_groups, dispersion)
    observed_u = randomized_nb_pit(
        counts,
        observed_fit.means,
        observed_fit.size,
        observed_rng.random(counts.size),
    )
    observed_stat = stratified_ks_stat(observed_u, groups, n_groups)

    fixed_stats = np.empty(boot, dtype=float)
    refitted_stats = np.empty(boot, dtype=float)
    start = time.perf_counter()
    for bootstrap_index in range(boot):
        bootstrap_counts = simulate_nb(rng, observed_fit.means, observed_fit.size)
        fixed_u = randomized_nb_pit(
            bootstrap_counts,
            observed_fit.means,
            observed_fit.size,
            rng.random(counts.size),
        )
        fixed_stats[bootstrap_index] = stratified_ks_stat(
            fixed_u, groups, n_groups
        )

        bootstrap_fit = fit_conditional_nb(
            bootstrap_counts, exposure, groups, n_groups, dispersion
        )
        refitted_u = randomized_nb_pit(
            bootstrap_counts,
            bootstrap_fit.means,
            bootstrap_fit.size,
            rng.random(counts.size),
        )
        refitted_stats[bootstrap_index] = stratified_ks_stat(
            refitted_u, groups, n_groups
        )

    fixed_exceedances = int(np.sum(fixed_stats >= observed_stat))
    refitted_exceedances = int(np.sum(refitted_stats >= observed_stat))
    lo, hi = monte_carlo_interval(refitted_exceedances, boot)
    result = {
        "gene": gene,
        "dispersion_model": dispersion,
        "cells": int(counts.size),
        "detection_rate": float(np.mean(counts > 0)),
        "fitted_size": float(np.median(observed_fit.size)),
        "fitted_size_min": float(np.min(observed_fit.size)),
        "fitted_size_max": float(np.max(observed_fit.size)),
        "mean_fitted_log_score": float(np.mean(nbinom.logpmf(
            counts, observed_fit.size, observed_fit.size / (observed_fit.size + observed_fit.means)
        ))),
        "observed_stratified_ks": observed_stat,
        "fixed_parameter_cutoff_95": float(
            np.quantile(fixed_stats, 0.95, method="higher")
        ),
        "refitted_cutoff_95": float(
            np.quantile(refitted_stats, 0.95, method="higher")
        ),
        "fixed_parameter_p": (1 + fixed_exceedances) / (boot + 1),
        "refitted_p": (1 + refitted_exceedances) / (boot + 1),
        "fixed_exceedances": fixed_exceedances,
        "refitted_exceedances": refitted_exceedances,
        "refitted_tail_mc_lower": lo,
        "refitted_tail_mc_upper": hi,
        "bootstrap_replicates": boot,
        "failed_fits": 0,
        "observed_seed": seed,
        "elapsed_seconds": float(time.perf_counter() - start),
    }
    if output_dir:
        destination = Path(output_dir) / "replicates"
        destination.mkdir(parents=True, exist_ok=True)
        stem = destination / f"{dispersion}_{gene}"
        np.savez_compressed(
            stem.with_suffix(".npz"), observed_statistic=observed_stat,
            fixed_statistics=fixed_stats, refitted_statistics=refitted_stats,
            observed_randomized_residuals=observed_u,
            fitted_size=observed_fit.size,
        )
        stem.with_suffix(".json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def monte_carlo_interval(k: int, b: int) -> tuple[float, float]:
    """Clopper-Pearson interval for the conditional exceedance probability."""
    return (
        float(beta.ppf(0.025, k, b - k + 1)) if k else 0.0,
        float(beta.ppf(0.975, k + 1, b - k)) if k < b else 1.0,
    )


def make_figure(
    qc: pd.DataFrame, diagnostics: pd.DataFrame, output: Path,
    comparison: pd.DataFrame | None = None,
) -> None:
    plt.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "font.family": "DejaVu Sans"})
    responder = "#287D8E"
    nonresponder = "#D97745"
    refitted = "#3A8062"
    fixed = "#B26155"
    ink = "#20252A"

    has_comparison = comparison is not None and "patient" in set(comparison.dispersion_model)
    figure = plt.figure(figsize=(7.0, 6.0), constrained_layout=True)
    grid = figure.add_gridspec(2, 2, height_ratios=(1, 1.15))
    axes = [figure.add_subplot(grid[0, 0] if has_comparison else grid[0, :]),
            figure.add_subplot(grid[0, 1] if has_comparison else grid[1, 0]),
            figure.add_subplot(grid[1, 0] if has_comparison else grid[1, 1])]
    if has_comparison:
        axes.append(figure.add_subplot(grid[1, 1]))

    colors = [
        responder if value == "Responder" else nonresponder
        for value in qc["response"]
    ]
    axes[0].bar(qc["sample"], qc["retained_cells"], color=colors, width=0.75)
    axes[0].set_ylabel("Cells retained after QC")
    axes[0].set_xlabel("Baseline PBMC sample")
    axes[0].tick_params(axis="x", rotation=45 if has_comparison else 0)
    axes[0].set_title("A  Retained cells", loc="left", fontsize=11, fontweight="bold", pad=9)
    axes[0].text(
        0.02,
        0.93,
        f"{qc['retained_cells'].sum():,} cells",
        transform=axes[0].transAxes,
        color=ink,
        fontsize=9,
    )
    axes[0].set_ylim(0, qc.retained_cells.max() * 1.25)
    axes[0].legend(
        handles=[
            Patch(facecolor=responder, label="Responder"),
            Patch(facecolor=nonresponder, label="Nonresponder"),
        ],
        frameon=False,
        fontsize=8.5,
        loc="upper left", bbox_to_anchor=(0, -0.29), ncol=2,
    )

    ordered = diagnostics.sort_values("observed_stratified_ks").reset_index(drop=True)
    y = np.arange(len(ordered))
    axes[1].hlines(
        y,
        ordered.refitted_cutoff_95,
        ordered.observed_stratified_ks,
        color="#C8CDD0",
        linewidth=1.2,
    )
    axes[1].scatter(
        ordered.observed_stratified_ks,
        y,
        marker="o",
        s=30,
        color=ink,
        label="Observed",
        zorder=3,
    )
    axes[1].scatter(
        ordered.fixed_parameter_cutoff_95,
        y,
        marker="^",
        s=34,
        color=fixed,
        label="Fixed 95%",
        zorder=3,
    )
    axes[1].scatter(
        ordered.refitted_cutoff_95,
        y,
        marker="s",
        s=28,
        color=refitted,
        label="Refitted 95%",
        zorder=3,
    )
    axes[1].set_yticks(y, ordered.gene)
    axes[1].set_xlabel(r"$\max_j\sqrt{n_j}D_j$")
    axes[1].legend(frameon=False, fontsize=8.5, loc="upper center", bbox_to_anchor=(.5, -0.34), ncol=2)
    axes[1].set_title("B  Shared-dispersion baseline", loc="left", fontsize=11, fontweight="bold", pad=9)

    positive_floor = 1.0 / (diagnostics.bootstrap_replicates.iloc[0] + 1)
    alternative = None
    if has_comparison:
        alternative = comparison[comparison.dispersion_model == "patient"].set_index("gene").loc[ordered.gene]
    first_p = ordered.refitted_p if has_comparison else ordered.fixed_parameter_p
    second_p = alternative.refitted_p if has_comparison else ordered.refitted_p
    axes[2].hlines(
        y,
        first_p,
        second_p,
        color="#C8CDD0",
        linewidth=1.2,
    )
    axes[2].scatter(
        first_p,
        y,
        marker="^",
        s=34,
        color=fixed,
        label="Shared dispersion" if has_comparison else "Fixed parameter",
        zorder=3,
    )
    axes[2].scatter(
        second_p,
        y,
        marker="s",
        s=28,
        color=refitted,
        label="Patient dispersion" if has_comparison else "Refitted",
        zorder=3,
    )
    axes[2].axvline(0.05, color="#6A6F73", linestyle="--", linewidth=1)
    axes[2].set_xscale("log")
    axes[2].set_xlim(positive_floor * 0.8, 1.05)
    axes[2].set_yticks(y, ordered.gene)
    axes[2].set_xlabel("Bootstrap p-value")
    if not has_comparison:
        axes[2].legend(frameon=False, fontsize=8.5, loc="upper center", bbox_to_anchor=(.5, -0.25), ncol=2)
    axes[2].set_title("C  Refitted tail probabilities", loc="left", fontsize=11, fontweight="bold", pad=9)

    if has_comparison:
        baseline_ratio = ordered.observed_stratified_ks / ordered.refitted_cutoff_95
        alternative_ratio = alternative.observed_stratified_ks / alternative.refitted_cutoff_95
        axes[3].hlines(y, baseline_ratio, alternative_ratio, color="#C8CDD0", linewidth=1.2)
        axes[3].scatter(baseline_ratio, y, marker="^", s=30, color=fixed, label="Shared dispersion", zorder=3)
        axes[3].scatter(alternative_ratio, y, marker="s", s=28, color=refitted, label="Patient dispersion", zorder=3)
        axes[3].set_yticks(y, ordered.gene)
        axes[3].axvline(1, color="#6A6F73", linestyle="--", linewidth=1)
        axes[3].set_xlim(0, max(baseline_ratio.max(), alternative_ratio.max()) * 1.08)
        axes[3].set_xlabel("Observed / refitted 95% cutoff")
        handles, labels = axes[3].get_legend_handles_labels()
        figure.legend(handles, labels, frameon=False, fontsize=9, loc="upper center", bbox_to_anchor=(.5, 0), ncol=2)
        axes[3].set_title("D  Relative discrepancy", loc="left", fontsize=11, fontweight="bold", pad=9)

    for axis in axes:
        axis.tick_params(labelsize=9.5)
        axis.xaxis.label.set_size(10)
        axis.yaxis.label.set_size(10)
        axis.spines[["top", "right"]].set_visible(False)
        axis.grid(axis="x", color="#E5E7E8", linewidth=0.7, zorder=0)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    figure.savefig(output.with_suffix(".png"), dpi=220, bbox_inches="tight")
    plt.close(figure)


def run(args: argparse.Namespace) -> None:
    data_dir = Path(args.data_dir)
    if args.download:
        download_missing(data_dir)

    samples: list[SampleData] = []
    qc_rows: list[dict[str, float | int | str]] = []
    retained_counts: dict[str, list[np.ndarray]] = {gene: [] for gene in PANEL}
    retained_library: list[np.ndarray] = []
    retained_groups: list[np.ndarray] = []

    for group, (accession, sample_name, response) in enumerate(SAMPLES):
        path = sample_path(data_dir, accession, sample_name)
        if not path.exists():
            raise FileNotFoundError(f"Missing {path}; rerun with --download")
        print(f"Reading {sample_name}", flush=True)
        sample = read_sample(path, sample_name, response)
        samples.append(sample)
        keep = qc_mask(sample)
        retained_library.append(sample.library_size[keep])
        retained_groups.append(np.full(np.sum(keep), group, dtype=np.int16))
        for gene in PANEL:
            retained_counts[gene].append(sample.panel_counts[gene][keep])

        qc_rows.append(
            {
                "accession": accession,
                "sample": sample_name,
                "response": response,
                "candidate_cells": int(sample.library_size.size),
                "retained_cells": int(np.sum(keep)),
                "retention_rate": float(np.mean(keep)),
                "median_umi_retained": float(np.median(sample.library_size[keep])),
                "median_genes_retained": float(np.median(sample.detected_genes[keep])),
                "median_mito_fraction_retained": float(
                    np.median(sample.mitochondrial_fraction[keep])
                ),
            }
        )

    qc = pd.DataFrame(qc_rows)
    library_size = np.concatenate(retained_library).astype(float)
    groups = np.concatenate(retained_groups)
    exposure = library_size / np.median(library_size)

    results_dir = Path(args.results_dir)
    results_dir.mkdir(parents=True, exist_ok=True)
    tasks = [dict(
        gene=gene, counts=np.concatenate(retained_counts[gene]).astype(np.int64),
        exposure=exposure, groups=groups, n_groups=len(SAMPLES), boot=args.boot,
        seed=args.seed + 1009 * gene_index, dispersion=model,
        output_dir=str(results_dir),
    ) for model in args.models for gene_index, gene in enumerate(PANEL)]
    diagnostics_rows = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futures = [pool.submit(analyze_gene, **task) for task in tasks]
        for future in as_completed(futures):
            result = future.result()
            diagnostics_rows.append(result)
            print(f"Completed {result['dispersion_model']} {result['gene']}: "
                  f"p={result['refitted_p']:.5g}, B={args.boot}", flush=True)

    diagnostics = pd.DataFrame(diagnostics_rows)
    diagnostics["panel_order"] = diagnostics.gene.map({g: i for i, g in enumerate(PANEL)})
    diagnostics = diagnostics.sort_values(["dispersion_model", "panel_order"]).drop(columns="panel_order")
    for model in args.models:
        mask = diagnostics.dispersion_model == model
        diagnostics.loc[mask, "holm_refitted_p"] = holm_adjust(
            diagnostics.loc[mask, "refitted_p"].to_numpy()
        )
    qc.to_csv(results_dir / "qc_summary.csv", index=False)
    diagnostics.to_csv(results_dir / "model_comparison.csv", index=False)
    baseline = diagnostics[diagnostics.dispersion_model == "shared"]
    if not baseline.empty:
        baseline.to_csv(results_dir / "refitted_margin_diagnostics.csv", index=False)
    source_paths = [sample_path(data_dir, item[0], item[1]) for item in SAMPLES]
    manifest = pd.DataFrame(
        {
            "accession": [item[0] for item in SAMPLES],
            "sample": [item[1] for item in SAMPLES],
            "url": [sample_url(item[0], item[1]) for item in SAMPLES],
            "bytes": [path.stat().st_size for path in source_paths],
            "sha256": [sha256sum(path) for path in source_paths],
        }
    )
    manifest.to_csv(results_dir / "source_manifest.csv", index=False)
    if not baseline.empty:
        make_figure(qc, baseline, Path(args.figure), diagnostics)

    print(qc.to_string(index=False))
    print(diagnostics.to_string(index=False))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-dir", default="data/genentech_imvigor210")
    parser.add_argument("--results-dir", default="results/genentech_imvigor210")
    parser.add_argument("--figure", default="figs/imvigor210_residual_diagnostics")
    parser.add_argument("--boot", type=int, default=199)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--models", nargs="+", choices=("shared", "patient"), default=["shared"])
    parser.add_argument("--seed", type=int, default=20260822)
    parser.add_argument("--download", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    run(parse_args())
