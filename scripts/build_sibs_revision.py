#!/usr/bin/env python3
"""Build the clean, highlighted, supplement, and response reading copies."""

from pathlib import Path
import argparse
import re
import shutil
import subprocess


ROOT = Path(__file__).resolve().parents[1]
BUILD = ROOT / "tmp/pdfs/revision_build"
PDF_OUT = ROOT / "output/pdf"
TEX_OUT = ROOT / "output/tex"
BASELINE = "fdcbb20"


def run(command: list[str], log_name: str) -> str:
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    (BUILD / log_name).write_text(result.stdout + result.stderr)
    if result.returncode:
        raise RuntimeError(f"{command[0]} failed; see {BUILD / log_name}")
    return result.stdout


def compile_tex(source: Path, output_directory: Path | None = None) -> Path:
    command = ["latexmk", "-pdf", "-interaction=nonstopmode", "-halt-on-error"]
    if output_directory:
        command.append(f"-outdir={output_directory}")
    command.append(str(source.relative_to(ROOT)))
    run(command, source.stem + "-build.txt")
    return (output_directory or source.parent) / (source.stem + ".pdf")


def body_tokens(source: str, marked: bool = False) -> str:
    source = re.sub(r"(?<!\\)%[^\n]*", "", source)
    body = source.split(r"\begin{document}", 1)[1]
    if marked:
        # Only unwrap annotations; keep every original TeX group intact.
        tokens = re.findall(r"\\[A-Za-z]+|\\.|[{}]|[^\\{}]+", body)
        stack = []
        accepted = []
        annotation_group = False
        for token in tokens:
            if token in (r"\DIFadd", r"\DIFaddFL"):
                annotation_group = True
            elif token in (r"\DIFaddbegin", r"\DIFaddend", r"\DIFaddbeginFL", r"\DIFaddendFL"):
                continue
            elif token.startswith(r"\DIFdel"):
                raise RuntimeError("Unexpected deletion in additions-only manuscript")
            elif token == "{":
                stack.append(annotation_group)
                if not annotation_group:
                    accepted.append(token)
                annotation_group = False
            elif token == "}":
                if not stack.pop():
                    accepted.append(token)
            else:
                accepted.append(token)
        if stack or annotation_group:
            raise RuntimeError("Unbalanced markup groups")
        body = "".join(accepted)
    return re.sub(r"\s+", "", body)


def expand_result_inputs(source: str) -> str:
    # Expand generated result fragments so their numbers receive revision marks.
    pattern = r"\\(?:input|inputresultrows)\{(results/(?:conditional_calibration|computational_reporting)/[A-Za-z_]+\.tex)\}"
    return re.sub(pattern, lambda match: (ROOT / match.group(1)).read_text(), source)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--public", action="store_true",
                        help="Build only the clean manuscript and supplement, without review files.")
    args = parser.parse_args()
    for directory in (BUILD, PDF_OUT, TEX_OUT):
        directory.mkdir(parents=True, exist_ok=True)

    source_names = {
        "ks_scdesign_sib": "sibs_major_revision_manuscript",
        "supplementary-appendix": "sibs_major_revision_supplement",
    }
    if not args.public:
        source_names["response_to_reviewers_sibs_major_revision"] = "sibs_major_revision_response_to_reviewers"
    for source_name, output_name in source_names.items():
        source = ROOT / (source_name + ".tex")
        compiled = compile_tex(source)
        shutil.copy2(compiled, PDF_OUT / (output_name + ".pdf"))
        shutil.copy2(source, TEX_OUT / (output_name + ".tex"))

    if args.public:
        print("Built the public manuscript and supplement PDFs; no review files required.")
        return

    archived_baseline = ROOT / "reproducibility/submitted_manuscript.tex"
    baseline_source = archived_baseline.read_text() if archived_baseline.exists() else run(
        ["git", "show", f"{BASELINE}:ks_scdesign_sib.tex"], "baseline-source.txt"
    )
    baseline_path = BUILD / "submitted_manuscript.tex"
    baseline_path.write_text(baseline_source)
    current_source = expand_result_inputs((ROOT / "ks_scdesign_sib.tex").read_text())
    current_path = BUILD / "current_manuscript_expanded.tex"
    current_path.write_text(current_source)
    marked = run(
        [
            "latexdiff", "--type=CFONT", "--no-del", "--math-markup=0",
            "--graphics-markup=none", "--no-label",
            str(baseline_path), str(current_path),
        ],
        "latexdiff.txt",
    )
    # Keep the original font metrics so highlighting does not widen tables.
    old_style = r"{\protect\color{blue} \sf #1}"
    if marked.count(old_style) != 1:
        raise RuntimeError("Unexpected latexdiff preamble; check highlight style")
    marked = marked.replace(old_style, r"{\protect\color{blue}#1}")
    marked_path = TEX_OUT / "sibs_major_revision_marked.tex"
    marked_path.write_text(marked)
    if body_tokens(marked, marked=True) != body_tokens(current_source):
        raise RuntimeError("Accepted marked text differs from clean manuscript")
    compiled = compile_tex(marked_path, BUILD)
    shutil.copy2(compiled, PDF_OUT / "sibs_major_revision_marked.pdf")
    print("Built four revision PDFs; marked text matches clean manuscript.")


if __name__ == "__main__":
    main()
