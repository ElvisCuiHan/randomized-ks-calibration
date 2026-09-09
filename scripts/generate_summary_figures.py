#!/usr/bin/env python3
"""Generate manuscript summary figures and result-based panels."""

from __future__ import annotations

import argparse
import csv
import os
import tempfile
from pathlib import Path

cache_dir = Path(
    os.environ.get("KS_FIGURE_MPL_CACHE", Path(tempfile.gettempdir()) / "ks_figure_mpl_cache")
)
cache_dir.mkdir(exist_ok=True)
os.environ.setdefault("MPLBACKEND", "Agg")
os.environ.setdefault("MPLCONFIGDIR", str(cache_dir))
os.environ.setdefault("XDG_CACHE_HOME", str(cache_dir))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch
from PIL import Image


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--results",
        default="results/ks_scdesign_extended_results.txt",
        help="Combined results file produced by the reproducibility scripts.",
    )
    parser.add_argument(
        "--figure-dir",
        default="figs",
        help="Directory where manuscript figures are written.",
    )
    return parser.parse_args()


def collect_block(lines: list[str], header: str) -> list[dict[str, str]]:
    for index, line in enumerate(lines):
        if line.strip() == header:
            rows: list[str] = []
            for row in lines[index + 1 :]:
                if not row.strip():
                    break
                rows.append(row)
            return list(csv.DictReader([header, *rows]))
    raise ValueError(f"Could not find results block header: {header}")


def parse_scdesign3(results_path: Path) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    lines = results_path.read_text(encoding="utf-8").splitlines()
    metadata = collect_block(
        lines,
        "dataset,cells,genes,cell_types,pseudotime_min,pseudotime_max",
    )
    summary = collect_block(
        lines,
        "dataset,fit,median_gene_ks_D,max_gene_ks_D,min_naive_uniform_p,"
        "max_abs_pairwise_residual_corr,max_pseudotime_quartile_ks",
    )
    return metadata, summary


def label_for_fit(fit: str) -> str:
    return {
        "moment_nb": "Moment NB",
        "scdesign3_conditional_nb": "scDesign3 NB",
    }.get(fit, fit)


def plot_scdesign3_results(
    metadata: list[dict[str, str]],
    summary: list[dict[str, str]],
    figure_dir: Path,
) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    dataset_order = [row["dataset"] for row in metadata]
    fit_order = ["moment_nb", "scdesign3_conditional_nb"]
    colors = {"moment_nb": "#B55239", "scdesign3_conditional_nb": "#006D77"}
    metrics = [
        ("median_gene_ks_D", "Median KS $D$", (0, 0.22)),
        ("max_abs_pairwise_residual_corr", "Max residual $|\\rho|$", (0, 1.0)),
        ("max_pseudotime_quartile_ks", "Max bin KS", (0, 0.24)),
    ]
    lookup = {(row["dataset"], row["fit"]): row for row in summary}
    meta_lookup = {row["dataset"]: row for row in metadata}
    dataset_labels = {
        "example_sce": "Human embryos",
        "example_count": "Mouse pancreas",
    }
    xlabels = [
        f"{dataset_labels.get(name, name.replace('_', ' '))}\n"
        f"{int(meta_lookup[name]['cells']):,} cells"
        for name in dataset_order
    ]

    fig, axes = plt.subplots(
        nrows=1,
        ncols=3,
        figsize=(6.7, 2.75),
        dpi=220,
        constrained_layout=True,
    )
    bar_width = 0.34
    x_positions = list(range(len(dataset_order)))
    for panel_index, (metric, title, ylim) in enumerate(metrics):
        ax = axes[panel_index]
        ax.text(
            -0.13,
            1.08,
            chr(ord("A") + panel_index),
            transform=ax.transAxes,
            fontsize=12,
            weight="bold",
            va="top",
        )
        for fit_index, fit in enumerate(fit_order):
            values = [
                float(lookup[(dataset, fit)][metric])
                for dataset in dataset_order
            ]
            offset = (fit_index - 0.5) * bar_width
            bars = ax.bar(
                [x + offset for x in x_positions],
                values,
                width=bar_width,
                color=colors[fit],
                alpha=0.9,
                label=label_for_fit(fit),
            )
            for bar, value in zip(bars, values):
                ax.text(
                    bar.get_x() + bar.get_width() / 2,
                    value + (ylim[1] * 0.025),
                    f"{value:.3f}",
                    ha="center",
                    va="bottom",
                    fontsize=7.5,
                )
        ax.set_title(title, loc="left", fontsize=9.0, weight="bold")
        ax.set_xticks(x_positions)
        ax.set_xticklabels(xlabels, fontsize=7.5)
        ax.tick_params(axis="y", labelsize=7.5)
        ax.set_ylim(*ylim)
        ax.set_ylabel("Residual discrepancy" if panel_index == 0 else "", fontsize=8)
        ax.spines[["top", "right"]].set_visible(False)
    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, frameon=False, fontsize=8, ncol=2,
               loc="upper center", bbox_to_anchor=(0.5, -0.015))
    fig.savefig(figure_dir / "scdesign3_residual_diagnostics.png", bbox_inches="tight")
    fig.savefig(figure_dir / "scdesign3_residual_diagnostics.pdf", bbox_inches="tight")
    plt.close(fig)


def plot_article_map(figure_dir: Path) -> None:
    figure_dir.mkdir(parents=True, exist_ok=True)
    source_path = figure_dir / "article_calibration_map_source.png"
    if not source_path.exists():
        raise FileNotFoundError(f"Missing source image for article map: {source_path}")

    image = Image.open(source_path)
    cropped = image.crop((0, 210, image.width, min(image.height, 850)))
    cropped.save(figure_dir / "article_calibration_map.png")
    cropped.convert("RGB").save(figure_dir / "article_calibration_map.pdf")


def plot_calibration_workflow(figure_dir: Path) -> None:
    """Draw Figure 1 at near-publication size with embedded vector text."""
    figure_dir.mkdir(parents=True, exist_ok=True)
    style = {
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
        "mathtext.fontset": "dejavusans",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    }
    with plt.rc_context(style):
        fig = plt.figure(figsize=(6.0, 5.2), dpi=180, facecolor="white")
        ax = fig.add_axes((0, 0, 1, 1))
        ax.set(xlim=(0, 6.0), ylim=(0, 5.2))
        ax.axis("off")
        ink, muted, connector = "#24313B", "#53616B", "#80909A"
        left, right = 0.46, 5.91
        content_width = right - left
        gap = 0.18
        node_width = (content_width - 2 * gap) / 3
        columns = [left + i * (node_width + gap) for i in range(3)]
        checked_labels = []

        def arrow(start, end, scale=8.5):
            ax.add_patch(FancyArrowPatch(
                start, end, arrowstyle="-|>", mutation_scale=scale,
                linewidth=0.85, color=connector, shrinkA=2, shrinkB=2,
                capstyle="round", joinstyle="round", zorder=1,
            ))

        def panel(x, y, width, height, face, accent):
            ax.add_patch(FancyBboxPatch(
                (x, y), width, height,
                boxstyle="round,pad=0,rounding_size=0.035",
                facecolor=face, edgecolor="#C9D3D9", linewidth=0.7,
                zorder=2,
            ))
            ax.plot([x + 0.055, x + width - 0.055],
                    [y + height, y + height], color=accent,
                    linewidth=1.65, solid_capstyle="round", zorder=3)
            return x, y, width, height

        def label(bounds, y, value, size=9.5, weight="normal", color=ink):
            x, _, width, _ = bounds
            artist = ax.text(
                x + width / 2, y, value, ha="center", va="center",
                fontsize=size, weight=weight, color=color,
                linespacing=1.25, zorder=4,
            )
            checked_labels.append((artist, bounds))

        stages = [
            (1, 4.98, "Specify the diagnostic target"),
            (5, 2.43, r"Refitted bootstrap: repeat $B$ times"),
            (6, 0.93, "Draw a component-specific conclusion"),
        ]
        for number, y, title in stages:
            ax.add_patch(Circle(
                (0.16, y), 0.115, facecolor="#F3F6F7",
                edgecolor="#BCC9D0", linewidth=0.8, zorder=3,
            ))
            ax.text(0.16, y, str(number), ha="center", va="center",
                    color=ink, fontsize=9.3, weight="bold", zorder=4)
            ax.text(left, y, title, va="center", fontsize=10.8,
                    weight="bold", color=ink)

        targets = [
            ("Marginal", "Gene-wise or\naggregated KS", "#EEF5FB", "#2D719F"),
            ("Covariate-indexed", "Residual variation over\npseudotime, batch, state", "#EEF7F3", "#287D6D"),
            ("Dependence", "Rank correlation or\nempirical copula", "#F5F1F8", "#795C96"),
        ]
        for x, (title, detail, face, accent) in zip(columns, targets):
            bounds = panel(x, 4.08, node_width, 0.68, face, accent)
            label(bounds, 4.56, title, size=10.0, weight="bold", color=accent)
            label(bounds, 4.28, detail, size=9.3)

        operations = [
            ("Fit the model", "$\\widehat\\theta=\\mathcal{A}(D)$\nConditional count model", "#287D6D", "#EEF7F3"),
            ("Randomized PIT", "$U_i=F_{\\widehat\\theta}(X_i^-\\mid Z_i)$\n$+V_i p_{\\widehat\\theta}(X_i\\mid Z_i)$", "#B66344", "#FCF3EE"),
            ("Compute statistic", "$S_{\\mathrm{obs}}=S(U,D)$\nTarget-specific statistic", "#795C96", "#F5F1F8"),
        ]
        for number, (x, (title, detail, accent, face)) in enumerate(zip(columns, operations), 2):
            ax.add_patch(Circle((x + 0.09, 3.84), 0.105, facecolor=accent, edgecolor="none"))
            ax.text(x + 0.09, 3.84, str(number), ha="center", va="center",
                    fontsize=9.0, color="white", weight="bold")
            ax.text(x + 0.25, 3.84, title, va="center", fontsize=9.4, weight="bold", color=accent)
            bounds = panel(x, 2.93, node_width, 0.68, face, accent)
            label(bounds, 3.27, detail, size=9.2)
        for first, second in zip(columns, columns[1:]):
            arrow((first + node_width, 3.30), (second, 3.30))
        ax.text(
            (left + right) / 2, 2.74,
            r"$V_i\overset{\mathrm{iid}}{\sim}\mathrm{Unif}(0,1)$, independent of the observed data",
            ha="center", va="center", fontsize=9.5, color=muted,
        )

        bootstrap = [
            ("Simulate", "$D_b^*\\sim P_{\\widehat\\theta}$\nPreserve the design"),
            ("Refit", "$\\widehat\\theta_b^*=\\mathcal{A}(D_b^*)$\nSame fitting algorithm"),
            ("Recompute", "Fresh independent PIT\nand the same statistic $S_b^*$"),
        ]
        for x, (title, detail) in zip(columns, bootstrap):
            bounds = panel(x, 1.42, node_width, 0.78, "#EFF6FA", "#2D719F")
            label(bounds, 2.01, title, size=10.0, weight="bold", color="#27658C")
            label(bounds, 1.68, detail, size=9.0)
        for first, second in zip(columns, columns[1:]):
            arrow((first + node_width, 1.81), (second, 1.81))
        ax.text((left + right) / 2, 1.20,
                r"Compare $S_{\mathrm{obs}}$ with the empirical distribution of $S_1^*,\ldots,S_B^*$",
                ha="center", va="center", fontsize=9.3, color=muted)
        ax.plot([left, right], [0.75, 0.75], color="#D4DDE2", linewidth=0.7)
        ax.text((left + right) / 2, 0.48,
                r"$\widehat p=\frac{1+\sum_{b=1}^{B}\mathbf{1}\{S_b^*\geq S_{\mathrm{obs}}\}}{B+1}$",
                ha="center", va="center", fontsize=12, color=ink)
        ax.text((left + right) / 2, 0.14,
                "Apply the prespecified multiplicity rule; report only the tested component.",
                ha="center", va="center", fontsize=9.0, color=muted)

        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        panel_extents = {}
        for artist, (x, y, width, height) in checked_labels:
            extent = artist.get_window_extent(renderer)
            low, high = ax.transData.transform(
                [(x + 0.04, y + 0.03), (x + width - 0.04, y + height - 0.03)]
            )
            if (extent.x0 < low[0] or extent.y0 < low[1]
                    or extent.x1 > high[0] or extent.y1 > high[1]):
                raise ValueError(f"Figure 1 label exceeds its panel: {artist.get_text()}")
            key = (x, y, width, height)
            for previous in panel_extents.get(key, []):
                if extent.overlaps(previous):
                    raise ValueError(f"Figure 1 labels overlap: {artist.get_text()}")
            panel_extents.setdefault(key, []).append(extent)

        for extension in ("pdf", "svg", "png"):
            fig.savefig(figure_dir / f"calibration_workflow_revised.{extension}",
                        dpi=450, bbox_inches="tight", pad_inches=0.035,
                        facecolor="white")
        plt.close(fig)


def main() -> None:
    args = parse_args()
    results_path = Path(args.results)
    figure_dir = Path(args.figure_dir)
    metadata, summary = parse_scdesign3(results_path)
    plot_article_map(figure_dir)
    plot_calibration_workflow(figure_dir)
    plot_scdesign3_results(metadata, summary, figure_dir)

    print("summary_figures")
    print("figure,path")
    print(f"article_calibration_map,{figure_dir / 'article_calibration_map.png'}")
    print(f"calibration_workflow_revised,{figure_dir / 'calibration_workflow_revised.png'}")
    print(f"scdesign3_residual_diagnostics,{figure_dir / 'scdesign3_residual_diagnostics.png'}")


if __name__ == "__main__":
    main()
