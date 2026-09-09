"""Refresh Wilson intervals from the stored synthetic experiment counts."""

import csv
import io
from pathlib import Path

from ks_scdesign_extended_simulations import rate_interval


ROOT = Path(__file__).resolve().parents[1]


def reconstruct_count(value, trials):
    printed = value.split(" ", 1)[0]
    count = round(float(printed) * trials)
    candidates = [k for k in range(max(0, count - 2), min(trials, count + 2) + 1)
                  if f"{k / trials:.3f}" == printed]
    if candidates != [count]:
        raise ValueError(f"Cannot recover a count unambiguously: {value}, M={trials}")
    return count


def main():
    path = ROOT / "results/ks_scdesign_extended_results.txt"
    blocks = path.read_text().strip().split("\n\n")
    settings = next(block for block in blocks if block.startswith("settings,"))
    arguments = dict(item.split("=", 1) for item in settings.split(",")[1:])
    mc, n = int(arguments["mc"]), int(arguments["n"])
    output, counts = [], []
    for block in blocks:
        name = block.splitlines()[0]
        if name == "conditional_residuals":
            output.append("conditional_calibration_source\n"
                          "Run scripts/ks_conditional_calibration.R; see results/conditional_calibration/summary.csv")
            continue
        if name not in {"fitted_nb_scale", "fitted_nb_margin", "copula_blind_spot"}:
            output.append(block)
            continue
        rows = list(csv.reader(block.splitlines()[1:]))
        header = rows[0]
        for row in rows[1:]:
            fields = range(4, 9) if name == "fitted_nb_scale" else range(1, len(row))
            for index in fields:
                if row[index] == "--":
                    continue
                trials = int(row[2] if index == 8 else row[1]) if name == "fitted_nb_scale" else mc
                count = reconstruct_count(row[index], trials)
                rate = count / trials
                lo, hi = rate_interval(rate, trials)
                row[index] = f"{rate:.3f} [{lo:.3f}; {hi:.3f}]"
                counts.append(dict(experiment=name, scenario=row[0],
                                   n=int(row[0]) if name == "fitted_nb_scale" else n * (2 if name == "copula_blind_spot" else 1),
                                   method=header[index].split("_with_95mc_ci")[0],
                                   trials=trials, rejections=count, rate=rate, lower=lo, upper=hi,
                                   interval="Wilson"))
        stream = io.StringIO()
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerows(rows)
        output.append(name + "\n" + stream.getvalue().rstrip())
    path.write_text("\n\n".join(output) + "\n")
    with (ROOT / "results/synthetic_rejection_counts.csv").open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(counts[0]))
        writer.writeheader()
        writer.writerows(counts)
    print(f"Refreshed {len(counts)} intervals; rejection counts preserved; legacy conditional block removed.")


if __name__ == "__main__":
    main()
