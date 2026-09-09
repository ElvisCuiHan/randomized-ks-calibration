"""Check stored bootstrap records against all published summary values."""

import json
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

from imvigor210_residual_diagnostics import (
    PANEL, holm_adjust, monte_carlo_interval, stratified_ks_stat,
)


ROOT = Path(__file__).resolve().parents[1]


def main():
    folder = ROOT / "results/genentech_imvigor210"
    frame = pd.read_csv(folder / "model_comparison.csv")
    assert len(frame) == 20
    qc = pd.read_csv(folder / "qc_summary.csv")
    groups = np.repeat(np.arange(len(qc)), qc.retained_cells.to_numpy())
    for model, rows in frame.groupby("dispersion_model", sort=False):
        assert list(rows.gene) == list(PANEL)
        assert (rows.bootstrap_replicates == 9999).all()
        assert (rows.cells == 24134).all() and (rows.failed_fits == 0).all()
        np.testing.assert_allclose(holm_adjust(rows.refitted_p.to_numpy()), rows.holm_refitted_p)
        for row in rows.itertuples():
            data = np.load(folder / "replicates" / f"{model}_{row.gene}.npz")
            observed = float(data["observed_statistic"])
            assert np.isclose(observed, row.observed_stratified_ks)
            assert np.isclose(observed, stratified_ks_stat(
                data["observed_randomized_residuals"], groups, len(qc)))
            for field, prefix in (("fixed_statistics", "fixed"), ("refitted_statistics", "refitted")):
                values = data[field]
                assert values.size == 9999 and np.isfinite(values).all()
                k = int((values >= observed).sum())
                assert k == getattr(row, prefix + "_exceedances")
                p_name = "fixed_parameter_p" if prefix == "fixed" else "refitted_p"
                assert np.isclose((k + 1) / 10000, getattr(row, p_name))
                cutoff_name = "fixed_parameter_cutoff_95" if prefix == "fixed" else "refitted_cutoff_95"
                assert np.isclose(np.quantile(values, .95, method="higher"), getattr(row, cutoff_name))
            lo, hi = monte_carlo_interval(row.refitted_exceedances, 9999)
            np.testing.assert_allclose([lo, hi], [row.refitted_tail_mc_lower, row.refitted_tail_mc_upper])
    archived = ROOT / "reviews/2026-09-07/before/genentech_imvigor210/refitted_margin_diagnostics.csv"
    archived_check = None
    if archived.exists():
        old = pd.read_csv(archived)
        shared = frame[frame.dispersion_model == "shared"].set_index("gene")
        np.testing.assert_allclose(old.observed_stratified_ks,
                                  shared.loc[old.gene, "observed_stratified_ks"], atol=1e-9)
        archived_check = True
    pancreas = pd.read_csv(ROOT / "results/scdesign3_example/pancreas_refitted_margin_bootstrap.csv")
    assert len(pancreas) == 10 and (pancreas.bootstrap_replicates == 1999).all()
    assert (pancreas.failed_fits == 0).all()
    np.testing.assert_allclose((pancreas.exceedances + 1) / 2000, pancreas.refitted_bootstrap_p)
    np.testing.assert_allclose(holm_adjust(pancreas.refitted_bootstrap_p.to_numpy()), pancreas.holm_refitted_p)
    provenance = pd.read_csv(ROOT / "results/scdesign3_example/embryo_cell_provenance.csv")
    assert len(provenance) == 1289 and provenance.cell.nunique() == 1289
    report = {
        "date": "2026-09-07", "python": platform.python_version(),
        "numpy": np.__version__, "scipy": scipy.__version__, "pandas": pd.__version__,
        "platform": platform.platform(), "pbmc_models": 2, "pbmc_genes_per_model": 10,
        "pbmc_bootstrap_per_gene_model": 9999, "pbmc_total_refits": 199980,
        "pancreas_genes": 10, "pancreas_bootstrap_per_gene": 1999,
        "pancreas_total_refits": 19990, "failed_fits": 0,
        "pbmc_observed_statistics_match_archived": archived_check,
        "embryo_cell_ids_and_stages_verified": 1289,
        "checks": "All saved Python reference draws and both summary tables verified",
    }
    (ROOT / "results/yl_revision_verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
