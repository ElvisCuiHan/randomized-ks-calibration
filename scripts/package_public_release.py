"""Export an allowlisted public snapshot without submission history or correspondence."""

import argparse
import hashlib
import json
import shutil
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile


ROOT = Path(__file__).resolve().parents[1]
VERSION = "1.0.0"
SCRIPT_NAMES = (
    "build_sibs_revision.py", "package_public_release.py", "generate_summary_figures.py",
    "ks_crossing_calculations.py", "ks_composite_simulation.py",
    "ks_scdesign_extended_simulations.py", "ks_scgtm_example.py",
    "ks_scdesign3_example.R", "imvigor210_residual_diagnostics.py", "session_info.R",
    "verify_scdesign3_provenance.R", "test_imvigor210_bootstrap.py",
    "test_yl_bootstrap.R", "verify_pancreas_results.R", "verify_revision_results.py",
    "check_revision_tables.py", "ks_conditional_calibration.R",
    "test_conditional_calibration.R", "refresh_simulation_intervals.py",
    "verify_conditional_calibration.R", "test_simulation_intervals.py",
    "write_conditional_tables.py", "report_computational_costs.py",
)
RESULT_DIRECTORIES = (
    "conditional_calibration", "computational_reporting", "scdesign3_example",
    "genentech_imvigor210", "scgtm_example",
)
RESULT_SUFFIXES = {".csv", ".json", ".txt", ".tex", ".rds", ".npz"}
FIGURE_NAMES = (
    "Zn.png", "calibration_workflow_revised.pdf", "calibration_workflow_revised.png",
    "calibration_workflow_revised.svg", "scgtm_residual_calibration.pdf",
    "scgtm_residual_calibration.png", "scgtm_calibration_workflow.pdf",
    "scgtm_calibration_workflow.png", "scdesign3_residual_diagnostics.pdf",
    "scdesign3_residual_diagnostics.png", "imvigor210_residual_diagnostics.pdf",
    "imvigor210_residual_diagnostics.png", "article_calibration_map_source.png",
    "article_calibration_map.pdf", "article_calibration_map.png",
)


def public_files():
    files = [ROOT / name for name in (
        "README.md", "LICENSE", "NOTICE.md", "CITATION.cff",
        "ks_scdesign_sib.tex", "supplementary-appendix.tex",
        "references.bib", "historical.bib", "ks_scdesign_sib.bbl",
        "supplementary-appendix.bbl", "output/pdf/sibs_major_revision_manuscript.pdf",
        "output/pdf/sibs_major_revision_supplement.pdf",
    )]
    files.extend(ROOT / "scripts" / name for name in SCRIPT_NAMES)
    files.extend(ROOT / "figs" / name for name in FIGURE_NAMES)
    files.extend(ROOT / "results" / name for name in (
        "ks_scdesign_extended_results.txt", "reproducibility_session_info.txt",
        "synthetic_rejection_counts.csv", "manuscript_table_verification.json",
        "yl_revision_verification.json",
    ))
    for directory in RESULT_DIRECTORIES:
        files.extend(p for p in (ROOT / "results" / directory).rglob("*")
                     if p.is_file() and p.suffix in RESULT_SUFFIXES
                     and not any(part.startswith(".") for part in p.relative_to(ROOT).parts))
    return sorted(set(files))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "output/public_release")
    args = parser.parse_args()
    destination = args.output_dir.resolve()
    if destination.exists():
        raise FileExistsError(f"Use a fresh export directory: {destination}")
    files = public_files()
    for path in files:
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"Missing or non-regular public input: {path}")
    destination.mkdir(parents=True)
    snapshot = destination / "repository"
    manifest = {}
    for path in files:
        relative = path.relative_to(ROOT)
        name = relative.as_posix()
        assert not any(part in {"reviews", "data", "tmp", ".git"} for part in relative.parts)
        assert "response_to_reviewers" not in name and "marked" not in name
        target = snapshot / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        manifest[name] = hashlib.sha256(target.read_bytes()).hexdigest()
    generated = {
        ".gitignore": "*.aux\n*.blg\n*.log\n*.out\n*.spl\n*.fdb_latexmk\n*.fls\n__pycache__/\n.DS_Store\ndata/\ntmp/\noutput/tex/\noutput/public_release*/\n",
        "reproducibility/release.json": json.dumps({
            "local_archive_version": VERSION, "date": "2026-09-09",
            "repository": "https://github.com/ElvisCuiHan/randomized-ks-calibration",
            "github_release_published": False,
            "scope": "Public manuscript, supplement, code and numerical evidence; no review correspondence.",
        }, indent=2) + "\n",
    }
    for name, contents in generated.items():
        target = snapshot / name
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(contents)
        manifest[name] = hashlib.sha256(target.read_bytes()).hexdigest()
    manifest_path = snapshot / "reproducibility/sha256.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    archive_path = destination / f"randomized-ks-calibration-v{VERSION}.zip"
    with ZipFile(archive_path, "w", ZIP_DEFLATED) as archive:
        for name in sorted(manifest):
            archive.write(snapshot / name, name)
        archive.write(manifest_path, "reproducibility/sha256.json")
    with ZipFile(archive_path) as archive:
        assert archive.testzip() is None
        assert set(archive.namelist()) == set(manifest) | {"reproducibility/sha256.json"}
        for name, expected in manifest.items():
            assert hashlib.sha256(archive.read(name)).hexdigest() == expected, name
    checksum = hashlib.sha256(archive_path.read_bytes()).hexdigest()
    (destination / "SHA256SUMS.txt").write_text(f"{checksum}  {archive_path.name}\n")
    print(f"Verified {len(manifest)} public files; {archive_path.stat().st_size:,} bytes")
    print(snapshot)
    print(archive_path)
    print(f"SHA-256: {checksum}")


if __name__ == "__main__":
    main()
