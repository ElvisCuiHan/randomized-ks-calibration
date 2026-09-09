"""Write manuscript fragments from the complete conditional calibration CSV."""

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
LABELS = {
    "linear_null": "Linear / linear",
    "periodic_misspecified": "Periodic / linear",
    "periodic_null": "Periodic / periodic",
    "spline_null": "Spline / spline",
}


def interval(row):
    return f"{row.rate:.3f} [{row.lower:.3f};{row.upper:.3f}]"


def main():
    folder = ROOT / "results/conditional_calibration"
    frame = pd.read_csv(folder / "summary.csv")
    assert len(frame) == 16 and set(frame.scenario) == set(LABELS)
    assert (frame.mc == 500).all() and (frame.boot == 199).all()
    assert (frame.failed_fits == 0).all()
    rows = []
    for scenario, label in LABELS.items():
        group = frame[frame.scenario == scenario].set_index(["statistic", "reference"])
        cells = [label, str(int(group.n.iloc[0]))]
        for key in (("marginal", "refitted"), ("conditional", "fixed"), ("conditional", "refitted")):
            cells.append(interval(group.loc[key]))
        rows.append(" & ".join(cells) + r"\\")
    (folder / "main_table_rows.tex").write_text("\n".join(rows) + "\n")
    rows = []
    for row in frame.itertuples():
        cells = [LABELS[row.scenario], row.statistic.capitalize(), row.reference.capitalize(),
                 str(row.n), f"{row.rejections}/{row.mc}", f"{row.rate:.3f}",
                 f"[{row.lower:.3f};{row.upper:.3f}]"]
        rows.append(" & ".join(cells) + r"\\")
    (folder / "supplement_table_rows.tex").write_text("\n".join(rows) + "\n")
    null = frame[(frame.scenario != "periodic_misspecified") & (frame.reference == "refitted")]
    marginal = null[null.statistic == "marginal"]
    conditional = null[null.statistic == "conditional"]
    alt = frame[(frame.scenario == "periodic_misspecified") & (frame.reference == "refitted")].set_index("statistic")
    text = (
        f"Across the three correctly specified cases, refitted marginal rejection\n"
        f"rates range from {100 * marginal.rate.min():.1f}\\% to {100 * marginal.rate.max():.1f}\\%, and\n"
        f"conditional rejection rates from {100 * conditional.rate.min():.1f}\\% to {100 * conditional.rate.max():.1f}\\%.\n"
        f"When the periodic component is omitted, conditional rejection is\n"
        f"{100 * alt.loc['conditional', 'rate']:.1f}\\%, compared with {100 * alt.loc['marginal', 'rate']:.1f}\\% for the marginal statistic.\n"
        "The two diagnostic targets thus have different sensitivity to the\n"
        "omitted trend in this setting. These are finite-sample\n"
        "checks of the specified regressions, with Monte Carlo uncertainty shown\n"
        "in the table, not a proof of validity for arbitrary fitted simulators.\n"
        "All 398,000 bootstrap refits and 2,000 observed fits succeeded.\n"
    )
    (folder / "main_results.tex").write_text(text)
    print("Generated conditional table fragments and result paragraph from all 16 summary rows.")


if __name__ == "__main__":
    main()
