#!/usr/bin/env python3
"""scGTM trend-fit residual diagnostics for the KS calibration manuscript.

This script uses the Python scGTM source package to fit pseudotime
trend models and then applies randomized PIT residual checks to the fitted
conditional count margins. The goal is to connect
the interpretable scGTM trend parameters to the residual calibration framework,
not to benchmark scGTM.

Example:
    SCGTM_SOURCE_DIR=/path/to/scGTM python3 scripts/ks_scgtm_example.py
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import logging
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import kstest, nbinom, poisson


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--scgtm-dir",
        default=os.environ.get("SCGTM_SOURCE_DIR"),
        help=(
            "Path to the scGTM source checkout. Either the repository root "
            "or its inner scGTM/ directory is accepted. Defaults to "
            "$SCGTM_SOURCE_DIR."
        ),
    )
    parser.add_argument(
        "--data",
        default=None,
        help=(
            "Optional CSV with a pseudotime column followed by gene count "
            "columns. If omitted, the scGTM demo CSV is used."
        ),
    )
    parser.add_argument(
        "--pseudotime-column",
        default="Time",
        help="Name of the pseudotime column in --data.",
    )
    parser.add_argument(
        "--genes",
        default="Gene1,Gene2,Gene3,Gene4,Gene5,Gene6,Gene7,Gene8,Gene9,Gene10",
        help="Comma-separated genes to fit.",
    )
    parser.add_argument(
        "--marginals",
        default="Poisson,NB",
        help="Comma-separated scGTM marginals to fit.",
    )
    parser.add_argument(
        "--iter",
        type=int,
        default=60,
        help="Number of PSO iterations passed to scGTM.",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=20260609,
        help="Seed for randomized PIT residuals.",
    )
    parser.add_argument(
        "--output-dir",
        default="results/scgtm_example",
        help="Directory where scGTM JSON and PNG outputs are written.",
    )
    parser.add_argument(
        "--figure-dir",
        default="figs",
        help="Directory where manuscript-ready scGTM summary figures are written.",
    )
    parser.add_argument(
        "--figure-genes",
        default="Gene5,Gene10",
        help="Comma-separated genes to show in the manuscript scGTM figure.",
    )
    parser.add_argument(
        "--skip-figures",
        action="store_true",
        help="Fit the models and write tables without regenerating manuscript figures.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Show scGTM optimization output.",
    )
    return parser.parse_args()


def resolve_scgtm_module_dir(scgtm_dir: str | None) -> Path:
    if not scgtm_dir:
        raise SystemExit(
            "scGTM source not found. Pass --scgtm-dir or set SCGTM_SOURCE_DIR."
        )

    base = Path(scgtm_dir).expanduser().resolve()
    candidates = [base, base / "scGTM"]
    for candidate in candidates:
        if (candidate / "scGTM.py").exists() and (candidate / "pseudotimeAPI.py").exists():
            return candidate
    raise SystemExit(f"Could not find scGTM.py under {base}")


def ks_stat_uniform(u: np.ndarray) -> float:
    values = np.sort(np.asarray(u, dtype=float))
    n = values.size
    upper = np.arange(1, n + 1, dtype=float) / n
    lower = np.arange(0, n, dtype=float) / n
    return float(max(np.max(upper - values), np.max(values - lower)))


def two_sample_ks(x: np.ndarray, y: np.ndarray) -> float:
    x = np.sort(np.asarray(x, dtype=float))
    y = np.sort(np.asarray(y, dtype=float))
    values = np.sort(np.concatenate([x, y]))
    fx = np.searchsorted(x, values, side="right") / x.size
    fy = np.searchsorted(y, values, side="right") / y.size
    return float(np.max(np.abs(fx - fy)))


def randomized_pit(
    counts: np.ndarray,
    means: np.ndarray,
    marginal: str,
    fit: dict[str, object],
    rng: np.random.Generator,
) -> np.ndarray:
    counts = np.asarray(counts, dtype=int)
    means = np.asarray(means, dtype=float)
    uniforms = rng.random(counts.size)

    if marginal == "Poisson":
        left = np.where(counts <= 0, 0.0, poisson.cdf(counts - 1, means))
        mass = poisson.pmf(counts, means)
    elif marginal == "NB":
        size = max(float(fit["phi"]), 1.0)
        prob = size / (size + means)
        left = np.where(counts <= 0, 0.0, nbinom.cdf(counts - 1, size, prob))
        mass = nbinom.pmf(counts, size, prob)
    else:
        raise ValueError("This residual script currently supports Poisson and NB.")

    return np.clip(left + uniforms * mass, 1e-12, 1.0 - 1e-12)


def max_quartile_ks(u: np.ndarray, pseudotime: np.ndarray) -> float:
    bins = pd.qcut(pseudotime, q=4, labels=False, duplicates="drop")
    max_stat = 0.0
    for i in sorted(set(bins)):
        for j in sorted(set(bins)):
            if j < i:
                max_stat = max(max_stat, two_sample_ks(u[bins == i], u[bins == j]))
    return max_stat


def fit_one(
    scgtm_module,
    link_fn,
    data: pd.DataFrame,
    pseudotime_column: str,
    gene: str,
    gene_index: int,
    marginal: str,
    iter_num: int,
    output_dir: Path,
    rng: np.random.Generator,
    verbose: bool,
) -> dict[str, object]:
    pseudotime = data[pseudotime_column].to_numpy(dtype=float)
    counts = np.floor(data[gene].to_numpy(dtype=float)).astype(int)

    call = lambda: scgtm_module.main(
        gene_index=gene_index,
        t=pseudotime,
        y1=counts,
        gene_name=gene,
        marginal=marginal,
        iter_num=iter_num,
        save_dir=str(output_dir) + "/",
        plot_args={
            "color": ["#005F73", "#0A9396", "#AE2012", "#9B2226"],
            "cmap": "viridis",
        },
    )
    if verbose:
        fit = call()
    else:
        with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
            fit = call()

    if int(fit.get("Transform", 0)) != 0:
        raise ValueError(
            "This residual diagnostic is defined on the original count scale; "
            "the selected scGTM fit used the package's valley transformation."
        )

    means = np.maximum(
        np.exp(link_fn(pseudotime, fit["mu"], fit["k1"], fit["k2"], fit["t0"])),
        0.1,
    )
    residuals = randomized_pit(counts, means, marginal, fit, rng)

    return {
        "gene": gene,
        "marginal": marginal,
        "cells": int(counts.size),
        "mean_count": float(np.mean(counts)),
        "zero_fraction": float(np.mean(counts == 0)),
        "negative_log_likelihood": float(fit["negative_log_likelihood"]),
        "mu": float(fit["mu"]),
        "k1": float(fit["k1"]),
        "k2": float(fit["k2"]),
        "t0": float(fit["t0"]),
        "phi": float(fit["phi"]) if marginal == "NB" else np.nan,
        "transform": int(fit.get("Transform", 0)),
        "t0_lower": float(fit["t0_lower"]),
        "t0_upper": float(fit["t0_upper"]),
        "ks_D": ks_stat_uniform(residuals),
        "ks_p_naive": float(kstest(residuals, "uniform").pvalue),
        "max_pseudotime_quartile_ks": max_quartile_ks(residuals, pseudotime),
        "_pseudotime": pseudotime,
        "_counts": counts,
        "_means": means,
        "_residuals": residuals,
    }


SUMMARY_COLUMNS = [
    "gene",
    "marginal",
    "cells",
    "mean_count",
    "zero_fraction",
    "negative_log_likelihood",
    "mu",
    "k1",
    "k2",
    "t0",
    "phi",
    "transform",
    "t0_lower",
    "t0_upper",
    "ks_D",
    "ks_p_naive",
    "max_pseudotime_quartile_ks",
]


def public_row(row: dict[str, object]) -> dict[str, object]:
    return {key: row[key] for key in SUMMARY_COLUMNS}


def json_ready(rows: list[dict[str, object]]) -> list[dict[str, object]]:
    ready_rows = []
    for row in rows:
        ready = {}
        for key, value in public_row(row).items():
            if isinstance(value, float) and np.isnan(value):
                ready[key] = None
            else:
                ready[key] = value
        ready_rows.append(ready)
    return ready_rows


def plot_scgtm_residual_calibration(
    rows: list[dict[str, object]],
    figure_dir: Path,
    figure_genes: list[str],
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure_dir.mkdir(parents=True, exist_ok=True)
    available = list(dict.fromkeys(str(row["gene"]) for row in rows))
    genes = [gene for gene in figure_genes if gene in available]
    if not genes:
        genes = available[:2]
    colors = {"Poisson": "#B55239", "NB": "#006D77"}
    styles = {"Poisson": "--", "NB": "-"}
    metric_labels = ["Uniform KS", "Quartile KS"]

    fig, axes = plt.subplots(
        nrows=len(genes),
        ncols=3,
        figsize=(11.4, 3.3 * len(genes)),
        dpi=220,
        constrained_layout=True,
    )
    if len(genes) == 1:
        axes = np.array([axes])

    for row_index, gene in enumerate(genes):
        gene_rows = [row for row in rows if row["gene"] == gene]
        base = gene_rows[0]
        pseudotime = np.asarray(base["_pseudotime"], dtype=float)
        counts = np.asarray(base["_counts"], dtype=float)
        order = np.argsort(pseudotime)

        ax = axes[row_index, 0]
        ax.text(
            -0.12,
            1.08,
            chr(ord("A") + row_index * 3),
            transform=ax.transAxes,
            fontsize=11,
            weight="bold",
            va="top",
        )
        ax.scatter(
            pseudotime,
            np.log1p(counts),
            s=12,
            c=pseudotime,
            cmap="viridis",
            alpha=0.72,
            linewidths=0,
        )
        for row in gene_rows:
            marginal = str(row["marginal"])
            means = np.asarray(row["_means"], dtype=float)
            ax.plot(
                pseudotime[order],
                np.log1p(means[order]),
                color=colors.get(marginal, "#333333"),
                linestyle=styles.get(marginal, "-"),
                linewidth=2.0,
                label=f"{marginal} fit",
            )
            ax.axvline(
                float(row["t0"]),
                color=colors.get(marginal, "#333333"),
                linestyle=styles.get(marginal, "-"),
                linewidth=1.2,
                alpha=0.75,
            )
        ax.set_title(f"{gene}: scGTM trend fit", loc="left", fontsize=10, weight="bold")
        ax.set_xlabel("Pseudotime")
        ax.set_ylabel("log(count + 1)")
        ax.legend(frameon=False, fontsize=8, loc="lower center")
        ax.spines[["top", "right"]].set_visible(False)

        ax = axes[row_index, 1]
        ax.text(
            -0.12,
            1.08,
            chr(ord("A") + row_index * 3 + 1),
            transform=ax.transAxes,
            fontsize=11,
            weight="bold",
            va="top",
        )
        grid = np.linspace(0, 1, 101)
        ax.plot(grid, grid, color="#777777", linewidth=1.0, linestyle=":", label="Uniform")
        for row in gene_rows:
            marginal = str(row["marginal"])
            residuals = np.sort(np.asarray(row["_residuals"], dtype=float))
            ecdf = np.arange(1, residuals.size + 1) / residuals.size
            ax.step(
                residuals,
                ecdf,
                where="post",
                color=colors.get(marginal, "#333333"),
                linewidth=1.8,
                label=f"{marginal}: D={float(row['ks_D']):.3f}",
            )
        ax.set_title("Randomized PIT calibration", loc="left", fontsize=10, weight="bold")
        ax.set_xlabel("Residual uniform scale")
        ax.set_ylabel("Empirical CDF")
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.legend(frameon=False, fontsize=8, loc="lower right")
        ax.spines[["top", "right"]].set_visible(False)

        ax = axes[row_index, 2]
        ax.text(
            -0.12,
            1.08,
            chr(ord("A") + row_index * 3 + 2),
            transform=ax.transAxes,
            fontsize=11,
            weight="bold",
            va="top",
        )
        x = np.arange(len(metric_labels))
        width = 0.32
        for offset, row in zip([-width / 2, width / 2], gene_rows):
            marginal = str(row["marginal"])
            values = [float(row["ks_D"]), float(row["max_pseudotime_quartile_ks"])]
            bars = ax.bar(
                x + offset,
                values,
                width=width,
                color=colors.get(marginal, "#333333"),
                alpha=0.88,
                label=marginal,
            )
            for bar, value in zip(bars, values):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    value + 0.012,
                    f"{value:.3f}",
                    ha="center",
                    va="bottom",
                    fontsize=7.5,
                )
        ax.set_title("Residual discrepancy", loc="left", fontsize=10, weight="bold")
        ax.set_xticks(x)
        ax.set_xticklabels(metric_labels)
        ax.set_ylim(0, 0.42)
        ax.set_ylabel("KS distance")
        ax.legend(frameon=False, fontsize=8, loc="upper right")
        ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle(
        "scGTM trend interpretability and randomized residual diagnostics",
        fontsize=12,
        weight="bold",
    )
    fig.savefig(figure_dir / "scgtm_residual_calibration.png", bbox_inches="tight")
    fig.savefig(figure_dir / "scgtm_residual_calibration.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_scgtm_workflow(figure_dir: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

    figure_dir.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10.4, 2.9), dpi=220)
    ax.set_axis_off()

    boxes = [
        ("Fit scGTM\ntrend margins", "pseudotime, counts\n$t_0,k_1,k_2,\\phi$"),
        ("Randomize\nPIT residuals", "$U_i=F(X_i-\\mid Z_i)+V_i p(X_i\\mid Z_i)$"),
        ("Use fitted-null\nreference", "refit or account for\nestimation term"),
        ("Report paired\ndiagnostics", "trend interpretation\n+ residual fit"),
    ]
    x_positions = [0.04, 0.295, 0.55, 0.805]
    colors = ["#E9D8A6", "#94D2BD", "#A9DEF9", "#F4A261"]

    for i, ((title, detail), x, color) in enumerate(zip(boxes, x_positions, colors)):
        box = FancyBboxPatch(
            (x, 0.32),
            0.17,
            0.38,
            boxstyle="round,pad=0.018,rounding_size=0.025",
            facecolor=color,
            edgecolor="#2F3E46",
            linewidth=1.1,
        )
        ax.add_patch(box)
        ax.text(x + 0.085, 0.575, title, ha="center", va="center", fontsize=10, weight="bold")
        ax.text(x + 0.085, 0.415, detail, ha="center", va="center", fontsize=8)
        if i < len(boxes) - 1:
            arrow = FancyArrowPatch(
                (x + 0.18, 0.51),
                (x_positions[i + 1] - 0.012, 0.51),
                arrowstyle="-|>",
                mutation_scale=14,
                linewidth=1.3,
                color="#2F3E46",
            )
            ax.add_patch(arrow)

    ax.text(
        0.5,
        0.12,
        "The residual diagnostic is descriptive unless the fitting step is included in the reference law.",
        ha="center",
        va="center",
        fontsize=9,
        color="#333333",
    )
    fig.savefig(figure_dir / "scgtm_calibration_workflow.png", bbox_inches="tight")
    fig.savefig(figure_dir / "scgtm_calibration_workflow.pdf", bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cache_dir = Path(
        os.environ.get("SCGTM_MPL_CACHE", Path(tempfile.gettempdir()) / "ks_scgtm_mpl_cache")
    )
    cache_dir.mkdir(exist_ok=True)
    os.environ.setdefault("MPLBACKEND", "Agg")
    os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
    os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))
    logging.getLogger("matplotlib").setLevel(logging.ERROR)
    logging.getLogger("pyswarms").setLevel(logging.ERROR)

    module_dir = resolve_scgtm_module_dir(args.scgtm_dir)
    sys.path.insert(0, str(module_dir))
    import scGTM as scgtm_module  # type: ignore
    from pseudotimeAPI import link as link_fn  # type: ignore

    data_path = Path(args.data) if args.data else module_dir / "Demo" / "simu_nb_scGTM_input.csv"
    data = pd.read_csv(data_path)
    display_data_path = data_path
    with contextlib.suppress(ValueError):
        display_data_path = data_path.resolve().relative_to(module_dir.parent.resolve())

    genes = [gene.strip() for gene in args.genes.split(",") if gene.strip()]
    figure_genes = [gene.strip() for gene in args.figure_genes.split(",") if gene.strip()]
    marginals = [m.strip() for m in args.marginals.split(",") if m.strip()]
    rng = np.random.default_rng(args.seed)
    np.random.seed(args.seed)

    rows = []
    for gene in genes:
        if gene not in data.columns:
            raise SystemExit(f"Gene column {gene!r} not found in {data_path}")
        gene_index = int(data.columns.get_loc(gene))
        for marginal in marginals:
            rows.append(
                fit_one(
                    scgtm_module=scgtm_module,
                    link_fn=link_fn,
                    data=data,
                    pseudotime_column=args.pseudotime_column,
                    gene=gene,
                    gene_index=gene_index,
                    marginal=marginal,
                    iter_num=args.iter,
                    output_dir=output_dir,
                    rng=rng,
                    verbose=args.verbose,
                )
            )

    if not args.skip_figures:
        figure_dir = Path(args.figure_dir)
        plot_scgtm_residual_calibration(rows, figure_dir, figure_genes)
        plot_scgtm_workflow(figure_dir)

    summary = pd.DataFrame([public_row(row) for row in rows])
    summary_path = output_dir / "scgtm_residual_summary.csv"
    summary.to_csv(summary_path, index=False)

    aggregate_rows = []
    if {"Poisson", "NB"}.issubset(set(summary["marginal"])):
        poisson = summary[summary["marginal"] == "Poisson"].set_index("gene")
        nb = summary[summary["marginal"] == "NB"].set_index("gene")
        common_genes = sorted(set(poisson.index).intersection(nb.index))
        if common_genes:
            ks_drop = poisson.loc[common_genes, "ks_D"] - nb.loc[common_genes, "ks_D"]
            bin_drop = (
                poisson.loc[common_genes, "max_pseudotime_quartile_ks"]
                - nb.loc[common_genes, "max_pseudotime_quartile_ks"]
            )
            nll_drop = (
                poisson.loc[common_genes, "negative_log_likelihood"]
                - nb.loc[common_genes, "negative_log_likelihood"]
            )
            aggregate_rows.append(
                {
                    "genes": len(common_genes),
                    "median_poisson_ks_D": float(poisson.loc[common_genes, "ks_D"].median()),
                    "median_nb_ks_D": float(nb.loc[common_genes, "ks_D"].median()),
                    "genes_with_lower_nb_ks_D": int((ks_drop > 0).sum()),
                    "median_ks_D_drop": float(ks_drop.median()),
                    "median_quartile_ks_drop": float(bin_drop.median()),
                    "median_nll_drop": float(nll_drop.median()),
                }
            )
    aggregate = pd.DataFrame(aggregate_rows)
    aggregate_path = output_dir / "scgtm_aggregate_summary.csv"
    aggregate.to_csv(aggregate_path, index=False)

    print("\nscgtm_example")
    print(
        "source,commit,source_dir,data,cells,genes,marginals,iter,seed,"
        "summary_csv,aggregate_csv"
    )
    commit = "unknown"
    head_file = module_dir.parent / ".git" / "HEAD"
    if head_file.exists():
        head = head_file.read_text(encoding="utf-8").strip()
        if head.startswith("ref: "):
            ref_file = module_dir.parent / ".git" / head.removeprefix("ref: ")
            if ref_file.exists():
                commit = ref_file.read_text(encoding="utf-8").strip()
        else:
            commit = head
    source_dir_display = "SCGTM_SOURCE_DIR" if os.environ.get("SCGTM_SOURCE_DIR") else str(module_dir.parent)
    print(
        f"ElvisCuiHan/scGTM,{commit},{source_dir_display},{display_data_path},"
        f"{len(data)},{len(genes)},{'|'.join(marginals)},{args.iter},"
        f"{args.seed},{summary_path},{aggregate_path}"
    )

    print("\nscgtm_reproduction_command")
    print(
        "SCGTM_SOURCE_DIR=/path/to/scGTM python3 scripts/ks_scgtm_example.py "
        f"--genes {','.join(genes)} --marginals {','.join(marginals)} "
        f"--iter {args.iter} --seed {args.seed}"
    )

    if not aggregate.empty:
        print("\nscgtm_aggregate_summary")
        print(
            "genes,median_poisson_ks_D,median_nb_ks_D,genes_with_lower_nb_ks_D,"
            "median_ks_D_drop,median_quartile_ks_drop,median_nll_drop,summary_csv"
        )
        agg = aggregate.iloc[0]
        print(
            f"{int(agg['genes'])},{agg['median_poisson_ks_D']:.3f},"
            f"{agg['median_nb_ks_D']:.3f},{int(agg['genes_with_lower_nb_ks_D'])},"
            f"{agg['median_ks_D_drop']:.3f},{agg['median_quartile_ks_drop']:.3f},"
            f"{agg['median_nll_drop']:.3f},{aggregate_path}"
        )

    print(
        "\ngene,marginal,cells,mean_count,zero_fraction,nll,mu,k1,k2,t0,"
        "phi,transform,t0_ci,ks_D,ks_p_naive,max_pseudotime_quartile_ks"
    )
    for row in rows:
        phi = "" if np.isnan(row["phi"]) else f"{row['phi']:.3f}"
        print(
            f"{row['gene']},{row['marginal']},{row['cells']},"
            f"{row['mean_count']:.3f},{row['zero_fraction']:.3f},"
            f"{row['negative_log_likelihood']:.3f},{row['mu']:.3f},"
            f"{row['k1']:.3f},{row['k2']:.3f},{row['t0']:.3f},"
            f"{phi},{row['transform']},[{row['t0_lower']:.3f};{row['t0_upper']:.3f}],"
            f"{row['ks_D']:.3f},{row['ks_p_naive']:.3g},"
            f"{row['max_pseudotime_quartile_ks']:.3f}"
        )

    json_path = output_dir / "scgtm_residual_summary.json"
    json_path.write_text(
        json.dumps(json_ready(rows), indent=2, allow_nan=False), encoding="utf-8"
    )


if __name__ == "__main__":
    main()
