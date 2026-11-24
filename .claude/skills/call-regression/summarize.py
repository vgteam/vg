#!/usr/bin/env python3
"""Print and save the accuracy summary of one scored run.

    summarize.py <run directory>

Reads what score.sh left in <run>/score/ and writes <run>/score/summary.json. Scores
from tools that did not run are left out.
"""

import csv
import json
import sys
from pathlib import Path


def aardvark(path: Path) -> dict:
    """Genotype-level F1 for all small variants, SNVs and indels."""
    names = {"ALL": "small", "Snv": "snv", "JointIndel": "indel"}
    out = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["comparison"] == "GT" and row["variant_type"] in names:
                key = names[row["variant_type"]]
                out[key] = {
                    "f1": float(row["metric_f1"]),
                    "truth_tp": int(row["truth_tp"]),
                    "truth_fn": int(row["truth_fn"]),
                    "query_fp": int(row["query_fp"]),
                }
    return out


def truvari(directory: Path) -> dict:
    """SV F1 from truvari bench, and after truvari refine where it ran."""
    out = {}
    bench = json.loads((directory / "summary.json").read_text())
    out["sv"] = {k: bench[k] for k in ("f1", "TP-base", "FN", "FP")}
    refined = directory / "refine.variant_summary.json"
    if refined.exists():
        out["sv_refined"] = {"f1": json.loads(refined.read_text())["f1"]}
    return out


def whatshap(path: Path) -> dict:
    """Phasing errors over all chromosomes compared.

    A flip is one site phased against both its neighbours, which whatshap's raw count
    reports as two switch errors. "switches" and "flips" are whatshap's decomposition
    of the raw count into the two.
    """
    raw = switches = flips = assessed = 0
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["chromosome"] == "ALL":
                continue
            assessed += int(row["all_assessed_pairs"])
            raw += int(row["all_switches"])
            s, f = row["all_switchflips"].split("/")
            switches += int(s)
            flips += int(f)
    return {"phasing": {"assessed_pairs": assessed, "switch_errors": raw,
                        "switches": switches, "flips": flips}}


def main() -> None:
    run = Path(sys.argv[1])
    score = run / "score"
    summary = {}
    summary.update(aardvark(score / "aardvark" / "summary.tsv"))
    if (score / "truvari" / "summary.json").exists():
        summary.update(truvari(score / "truvari"))
    if (score / "whatshap.tsv").exists():
        summary.update(whatshap(score / "whatshap.tsv"))
    (score / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    cells = [f"{k} F1 {v['f1']:.6f}" for k, v in summary.items() if "f1" in v]
    if "phasing" in summary:
        p = summary["phasing"]
        cells.append(f"switches {p['switches']} flips {p['flips']} of {p['assessed_pairs']} pairs")
    print(f"{run}: " + ", ".join(cells))


if __name__ == "__main__":
    main()
