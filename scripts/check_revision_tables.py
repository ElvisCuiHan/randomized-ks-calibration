"""Check gene-level, simulation and computational tables against stored results."""

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from ks_scdesign_extended_simulations import rate_interval
from write_conditional_tables import LABELS
from build_sibs_revision import expand_result_inputs


ROOT = Path(__file__).resolve().parents[1]


def table_rows(source, anchor):
    text = source[source.index(anchor):]
    text = text[text.index(r"\begin{tabular}"):text.index(r"\end{tabular}")]
    rows = {}
    for line in text.splitlines():
        if "&" not in line or not line.rstrip().endswith(r"\\"):
            continue
        fields = [field.strip() for field in line.rstrip()[:-2].split("&")]
        gene = re.sub(r"\\textit\{([^}]+)\}", r"\1", fields[0])
        if gene == "Gene":
            continue
        rows[gene] = fields[1:]
    return rows


def check_fields(fields, expected, decimals):
    assert len(fields) == len(expected), (fields, expected)
    for field, value, places in zip(fields, expected, decimals):
        assert abs(float(field.strip("$")) - value) <= 0.5 * 10 ** (-places) + 1e-10, (field, value)


def interval_numbers(field):
    values = [float(value) for value in re.findall(r"\d+\.\d+", field)]
    assert len(values) == 3, field
    return values


def check_simulation_tables(main_tex, supplement):
    counts = pd.read_csv(ROOT / "results/synthetic_rejection_counts.csv")
    assert len(counts) == 27
    for row in counts.itertuples():
        assert np.isclose(row.rate, row.rejections / row.trials, rtol=1e-12, atol=1e-15)
        np.testing.assert_allclose([row.lower, row.upper], rate_interval(row.rate, row.trials), atol=1e-12)
        assert row.upper > row.lower
    scale = table_rows(main_tex, r"\label{tab:nb-scale}")
    methods = ["known_margin", "naive_ks", "fixed_parameter_bootstrap",
               "refitted_bootstrap", "zero_inflated_power"]
    for n in (250, 1000, 10000):
        selected = counts[(counts.experiment == "fitted_nb_scale") & (counts.n == n)].set_index("method")
        check_fields(scale[f"{n:,}"][:5], [selected.loc[method, "rate"] for method in methods], [3] * 5)
    copula = table_rows(main_tex, r"\label{tab:scdesign-simulations}")
    names = {
        r"Same ($\rho=0.65$)": "same copula (rho=0.65)",
        r"Wrong ($\rho=0.50$)": "moderate wrong copula (rho=0.50)",
        r"Independent ($\rho=0$)": "rho=0 independence copula",
    }
    methods = ["gene1_marginal_ks", "gene2_marginal_ks", "residual_correlation"]
    for label, scenario in names.items():
        selected = counts[(counts.experiment == "copula_blind_spot") & (counts.scenario == scenario)].set_index("method")
        for cell, method in zip(copula[label], methods):
            row = selected.loc[method]
            np.testing.assert_allclose(interval_numbers(cell), [row.rate, row.lower, row.upper], atol=.00050001, rtol=0)
    folder = ROOT / "results/conditional_calibration"
    frame = pd.read_csv(folder / "summary.csv")
    assert len(frame) == 16
    main_rows = (folder / "main_table_rows.tex").read_text()
    assert r"\inputresultrows{results/conditional_calibration/main_table_rows.tex}" in main_tex
    assert r"\inputresultrows{results/conditional_calibration/supplement_table_rows.tex}" in supplement
    for line, (scenario, label) in zip(main_rows.splitlines(), LABELS.items()):
        fields = [item.strip() for item in line[:-2].split("&")]
        assert fields[0] == label and len(fields) == 5
        selected = frame[frame.scenario == scenario].set_index(["statistic", "reference"])
        assert int(fields[1]) == selected.n.iloc[0]
        keys = [("marginal", "refitted"), ("conditional", "fixed"), ("conditional", "refitted")]
        for cell, key in zip(fields[2:], keys):
            row = selected.loc[key]
            np.testing.assert_allclose(interval_numbers(cell), [row.rate, row.lower, row.upper], atol=.00050001, rtol=0)
    lines = (folder / "supplement_table_rows.tex").read_text().splitlines()
    assert len(lines) == 16 and len(main_rows.splitlines()) == 4
    for line, row in zip(lines, frame.itertuples()):
        fields = [item.strip() for item in line[:-2].split("&")]
        assert fields[:5] == [LABELS[row.scenario], row.statistic.capitalize(), row.reference.capitalize(),
                               str(row.n), f"{row.rejections}/{row.mc}"]
        np.testing.assert_allclose(interval_numbers(fields[5] + " " + fields[6]),
                                   [row.rate, row.lower, row.upper], atol=.00050001, rtol=0)


def check_computational_tables(main_tex, supplement):
    folder = ROOT / "results/computational_reporting"
    repeats = pd.read_csv(folder / "lightweight_repetitions.csv")
    summary = pd.read_csv(folder / "lightweight_summary.csv")
    assert len(repeats) == 20 and len(summary) == 4
    rows = table_rows(expand_result_inputs(main_tex), r"\label{tab:scalability}")
    for row in summary.itertuples():
        group = repeats[repeats.n == row.n]
        assert len(group) == 5 and set(group["repeat"]) == set(range(1, 6))
        obs, refit = group.observed_seconds.median(), group.one_refit_seconds.median()
        projection, memory = (obs + 199 * refit) / 60, 24 * row.n / 1024 ** 2
        np.testing.assert_allclose([row.observed_seconds, row.one_refit_seconds,
                                   row.projected_199_serial_minutes, row.streaming_working_arrays_mb],
                                  [obs, refit, projection, memory])
        check_fields(rows[f"{row.n:,}"], [10, obs, refit, projection, memory], [0, 3, 3, 3, 2])
    frame = pd.read_csv(folder / "task_runtime_summary.csv")
    assert len(frame) == 7 and frame.bootstrap_refits.sum() == 617970
    source = pd.read_csv(folder / "conditional_task_times.csv")
    assert len(source) == 2000 and not source.duplicated(["scenario", "task"]).any()
    sources = {label: source[source.scenario == key].elapsed_seconds for key, label in LABELS.items()}
    sources["Pancreas"] = pd.read_csv(ROOT / "results/scdesign3_example/pancreas_refitted_margin_bootstrap.csv").elapsed_seconds
    pbmc = pd.read_csv(ROOT / "results/genentech_imvigor210/model_comparison.csv")
    for model in ("shared", "patient"):
        sources[f"PBMC / {model}"] = pbmc[pbmc.dispersion_model == model].elapsed_seconds
    rows = table_rows(expand_result_inputs(supplement), r"\label{tab:computational-workloads}")
    for row in frame.itertuples():
        times = sources[row.workload]
        assert row.tasks == len(times) and row.bootstrap_refits == row.tasks * row.boot
        np.testing.assert_allclose([row.median_task_seconds, row.summed_task_seconds],
                                  [times.median(), times.sum()])
        fields = [field.replace(",", "") for field in rows[row.workload]]
        check_fields(fields, [row.cells, row.tasks, row.boot, row.bootstrap_refits,
                              times.median(), times.sum() / 60], [0, 0, 0, 0, 2, 2])


def main():
    main_tex = (ROOT / "ks_scdesign_sib.tex").read_text()
    supplement = (ROOT / "supplementary-appendix.tex").read_text()
    pancreas = pd.read_csv(ROOT / "results/scdesign3_example/pancreas_refitted_margin_bootstrap.csv")
    pbmc = pd.read_csv(ROOT / "results/genentech_imvigor210/model_comparison.csv")
    baseline = pbmc[pbmc.dispersion_model == "shared"].set_index("gene")
    extension = pbmc[pbmc.dispersion_model == "patient"].set_index("gene")
    rows = table_rows(main_tex, r"\label{tab:pancreas-refit}")
    assert set(rows) == set(pancreas.gene)
    for row in pancreas.itertuples():
        check_fields(rows[row.gene], [row.observed_ks_D, row.naive_uniform_p,
                     row.refitted_bootstrap_p, row.holm_refitted_p, row.bootstrap_q95], [4] * 5)
    rows = table_rows(supplement, "Shared-dispersion IMvigor210")
    assert set(rows) == set(baseline.index)
    for gene, row in baseline.iterrows():
        check_fields(rows[gene], [row.detection_rate, row.fitted_size,
                     row.observed_stratified_ks, row.fixed_parameter_cutoff_95,
                     row.refitted_cutoff_95, row.fixed_parameter_p, row.refitted_p,
                     row.holm_refitted_p], [3] * 5 + [4] * 3)
    rows = table_rows(supplement, "Exploratory patient-specific-dispersion comparison.")
    assert set(rows) == set(extension.index)
    for gene, row in extension.iterrows():
        old = baseline.loc[gene]
        check_fields(rows[gene], [old.observed_stratified_ks, row.observed_stratified_ks,
                     row.refitted_cutoff_95, row.refitted_p, row.holm_refitted_p,
                     row.mean_fitted_log_score - old.mean_fitted_log_score], [3, 3, 3, 4, 4, 6])
    rows = table_rows(supplement, "Conditional Monte Carlo uncertainty for the IMvigor210")
    for gene in baseline.index:
        fields = rows[gene]
        for i, frame in enumerate((baseline, extension)):
            row = frame.loc[gene]
            assert int(fields[2 * i]) == row.refitted_exceedances
            interval = [float(x) for x in fields[2 * i + 1].strip("[]").split(",")]
            np.testing.assert_allclose(interval, [row.refitted_tail_mc_lower,
                                       row.refitted_tail_mc_upper], atol=5.01e-7, rtol=0)
    rows = table_rows(supplement, "Pancreatic bootstrap exceedances and conditional")
    for row in pancreas.itertuples():
        fields = rows[row.gene]
        check_fields(fields[:2], [row.exceedances, row.refitted_bootstrap_p], [0, 4])
        interval = [float(x) for x in fields[2].strip("[]").split(",")]
        np.testing.assert_allclose(interval, [row.tail_mc_lower, row.tail_mc_upper], atol=5.01e-7, rtol=0)
    check_simulation_tables(main_tex, supplement)
    check_computational_tables(main_tex, supplement)
    report = {"date": "2026-09-09", "tables_checked": 10, "gene_rows_checked": 50,
              "simulation_rows_checked": 26, "synthetic_intervals_checked": 27,
              "runtime_rows_checked": 11, "lightweight_repetitions_checked": 20,
              "status": "passed", "rounding": "Verified against each displayed decimal precision"}
    (ROOT / "results/manuscript_table_verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
