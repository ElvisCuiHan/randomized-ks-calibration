"""Summarize saved task timers and benchmark the existing lightweight workflow."""

import argparse
from dataclasses import asdict
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import subprocess
import time

import numpy as np
import pandas as pd
import scipy

from ks_scdesign_extended_simulations import benchmark_streaming_scalability


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results/computational_reporting"
SCENARIOS = {
    "linear_null": "Linear / linear",
    "periodic_misspecified": "Periodic / linear",
    "periodic_null": "Periodic / periodic",
    "spline_null": "Spline / spline",
}


def summarize_tasks():
    extraction = r'''
    paths <- list.files("results/conditional_calibration", pattern="^trial-.*[.]rds$",
                        recursive=TRUE, full.names=TRUE)
    rows <- lapply(paths, function(path) {
      x <- readRDS(path)
      data.frame(scenario=x$scenario, task=x$trial, n=x$n, boot=x$boot,
                 elapsed_seconds=x$elapsed_seconds, failed_fits=x$failed_fits)
    })
    write.csv(do.call(rbind, rows), stdout(), row.names=FALSE)
    '''
    output = subprocess.check_output(["Rscript", "-e", extraction], cwd=ROOT, text=True)
    conditional = pd.read_csv(io.StringIO(output))
    assert len(conditional) == 2000 and not conditional.duplicated(["scenario", "task"]).any()
    assert (conditional.failed_fits == 0).all() and (conditional.boot == 199).all()
    conditional.to_csv(OUT / "conditional_task_times.csv", index=False)
    rows = []

    def append(label, cells, tasks, boot, times, scope):
        assert len(times) == tasks and np.isfinite(times).all() and (times > 0).all()
        rows.append(dict(workload=label, cells=int(cells), tasks=int(tasks), boot=int(boot),
                         bootstrap_refits=int(tasks * boot), median_task_seconds=float(np.median(times)),
                         min_task_seconds=float(np.min(times)), max_task_seconds=float(np.max(times)),
                         summed_task_seconds=float(np.sum(times)), timing_scope=scope))

    for key, label in SCENARIOS.items():
        group = conditional[conditional.scenario == key]
        assert len(group) == 500
        append(label, group.n.iloc[0], 500, 199, group.elapsed_seconds,
               "One outer sample: design, observed fit and bootstrap; excludes record serialization")
    pancreas = pd.read_csv(ROOT / "results/scdesign3_example/pancreas_refitted_margin_bootstrap.csv")
    assert len(pancreas) == 10 and (pancreas.bootstrap_replicates == 1999).all()
    assert (pancreas.failed_fits == 0).all()
    append("Pancreas", pancreas.cells.iloc[0], 10, 1999, pancreas.elapsed_seconds,
           "One gene: observed fit, bootstrap, RDS output and tail summaries; excludes data loading")
    pbmc = pd.read_csv(ROOT / "results/genentech_imvigor210/model_comparison.csv")
    for model in ("shared", "patient"):
        group = pbmc[pbmc.dispersion_model == model]
        assert len(group) == 10 and (group.bootstrap_replicates == 9999).all()
        assert (group.failed_fits == 0).all()
        append(f"PBMC / {model}", group.cells.iloc[0], 10, 9999, group.elapsed_seconds,
               "One gene-model: bootstrap and summaries; excludes observed fit, loading, QC and output")
    frame = pd.DataFrame(rows)
    assert frame.bootstrap_refits.sum() == 617970
    frame.to_csv(OUT / "task_runtime_summary.csv", index=False)
    lines = []
    for row in frame.itertuples():
        fields = [row.workload, f"{row.cells:,}", f"{row.tasks:,}", f"{row.boot:,}",
                  f"{row.bootstrap_refits:,}", f"{row.median_task_seconds:.2f}",
                  f"{row.summed_task_seconds / 60:.2f}"]
        lines.append(" & ".join(fields) + r"\\")
    (OUT / "task_table_rows.tex").write_text("\n".join(lines) + "\n")
    return frame


def run_lightweight_benchmark():
    for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS"):
        if os.environ.get(variable) != "1":
            raise RuntimeError(f"Set {variable}=1 before starting Python")
    base_seed = 20260909
    benchmark_streaming_scalability(np.random.default_rng(base_seed), [10000], 10)
    rows = []
    start = time.perf_counter()
    for repeat in range(1, 6):
        rng = np.random.default_rng(base_seed + repeat)
        for result in benchmark_streaming_scalability(rng, [1000, 10000, 100000, 1000000], 10):
            rows.append(dict(repeat=repeat, seed=base_seed + repeat, **asdict(result)))
        print(f"Lightweight benchmark repetition {repeat}/5 complete", flush=True)
    frame = pd.DataFrame(rows)
    frame.to_csv(OUT / "lightweight_repetitions.csv", index=False)
    metadata = dict(date_utc=datetime.now(timezone.utc).isoformat(), repeats=5, seed=base_seed,
                    sizes=[1000, 10000, 100000, 1000000], genes=10,
                    warmup="One n=10000, G=10 pass excluded from reported repetitions",
                    python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                    pandas=pd.__version__, platform=platform.platform(), workers=1,
                    environment={key: os.environ[key] for key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS")},
                    wall_seconds_excluding_imports_and_warmup=time.perf_counter() - start,
                    memory="24*n bytes analytic three-array baseline; peak RSS not measured",
                    source_sha256=hashlib.sha256((ROOT / "scripts/ks_scdesign_extended_simulations.py").read_bytes()).hexdigest(),
                    hardware=json.loads((OUT / "hardware_inventory.json").read_text()))
    (OUT / "lightweight_environment.json").write_text(json.dumps(metadata, indent=2) + "\n")


def write_lightweight_table():
    frame = pd.read_csv(OUT / "lightweight_repetitions.csv")
    assert len(frame) == 20 and set(frame["repeat"]) == set(range(1, 6))
    rows, lines = [], []
    for n, group in frame.groupby("n", sort=True):
        assert len(group) == 5 and (group.genes == 10).all()
        observed = float(group.observed_seconds.median())
        refit = float(group.one_refit_seconds.median())
        projected = (observed + 199 * refit) / 60
        memory = 24 * n / 1024 ** 2
        rows.append(dict(n=n, genes=10, observed_seconds=observed, one_refit_seconds=refit,
                         projected_199_serial_minutes=projected, streaming_working_arrays_mb=memory,
                         observed_min=group.observed_seconds.min(), observed_max=group.observed_seconds.max(),
                         refit_min=group.one_refit_seconds.min(), refit_max=group.one_refit_seconds.max()))
        fields = [f"{n:,}", "10", f"{observed:.3f}", f"{refit:.3f}", f"{projected:.3f}", f"{memory:.2f}"]
        lines.append(" & ".join(fields) + r"\\")
    pd.DataFrame(rows).to_csv(OUT / "lightweight_summary.csv", index=False)
    (OUT / "lightweight_table_rows.tex").write_text("\n".join(lines) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark", action="store_true", help="Run and overwrite only the lightweight timing experiment")
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    tasks = summarize_tasks()
    if args.benchmark:
        run_lightweight_benchmark()
    write_lightweight_table()
    print(tasks.to_string(index=False))


if __name__ == "__main__":
    main()
