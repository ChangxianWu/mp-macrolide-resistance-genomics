#!/usr/bin/env python3
"""Run the analysis pipeline stage by stage.

The repository keeps all pipeline scripts flat in ``scripts/``; this runner is
the single entry point and documents the order. Every script resolves the
repository root as the parent of ``scripts/``.

Examples
--------
    python run_pipeline.py --list
    python run_pipeline.py --dry-run
    python run_pipeline.py --stage 04_typing
    python run_pipeline.py --from 03_reference --to 05_coresnp_tree
"""
from __future__ import annotations

import argparse
import os
import runpy
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCRIPTS = ROOT / "scripts"

# Stage name -> [(script, extra arguments), ...] in execution order.
STAGES: list[tuple[str, list[tuple[str, list[str]]]]] = [
    ("01_panel_metadata", [
        ("01_build_manifest.py", []),
        ("01b_ena_manifest.py", []),
    ]),
    ("02_acquisition", [
        ("02_download_fasta.py", ["--included-only", "--workers", "6"]),
    ]),
    ("03_reference", [
        ("03e_calibrate_sites.py", []),
    ]),
    ("04_typing", [
        ("04_type_23S.py", []),
        ("06_summary_tables.py", []),
        ("07_w2qc.py", []),
        ("08_p1_type.py", []),
        ("09_mlst.py", []),
        ("10_matrix.py", []),
    ]),
    ("05_coresnp_tree", [
        ("11_p1_cluster.py", []),
        ("12_mash_nj.py", []),
        ("14_coresnp.py", []),
        ("15_snp_matrix.py", []),
    ]),
    ("06_temporal", [
        ("15b_clock.py", []),
        ("16_export_aln.py", []),
        ("17_datetree.py", []),
        ("18_export_full.py", []),
        ("19_treetime.py", []),
        ("19b_treetime_report.py", []),
        ("20_resistance_origin.py", []),
        ("21_epi_dynamics.py", []),
    ]),
    ("07_group_structure", [
        ("23_epidemic_clade.py", []),
        ("26_fulltree_ancestry.py", []),
        ("32_tree_perturbation_fitch.py", []),
        ("33_cluster_substructure.py", []),
    ]),
    ("08_normalization", [
        ("28_country_norm.py", []),
    ]),
    ("09_figures_tables", [
        ("22_figures_en.py", []),
        ("24_manuscript_tables.py", []),
        ("25_supplementary.py", []),
        ("29_fig3_redesign.py", []),
        ("30_figS1_threshold.py", []),
        ("31_figS2_composition.py", []),
        ("34_fig5_temporal_series.py", []),
    ]),
]

NAMES = [name for name, _ in STAGES]


def _index(name: str) -> int:
    if name in NAMES:
        return NAMES.index(name)
    if name.isdigit():
        n = int(name)
        if 1 <= n <= len(NAMES):
            return n - 1
    raise SystemExit(f"Unknown stage {name!r}. Available: {', '.join(NAMES)}")


def _select(start: str | None, stop: str | None, only: str | None):
    if only:
        i = _index(only)
        return STAGES[i:i + 1]
    lo = _index(start) if start else 0
    hi = _index(stop) if stop else len(STAGES) - 1
    if lo > hi:
        raise SystemExit("--from must not be after --to")
    return STAGES[lo:hi + 1]


def _print_plan(stages) -> None:
    for name, steps in stages:
        print(f"[{name}]")
        for script, args in steps:
            path = SCRIPTS / script
            flag = "" if path.is_file() else "   (missing!)"
            print("  python " + " ".join([f"scripts/{script}", *args]) + flag)
    print("\nDry run: nothing executed.")


def _run(stages) -> int:
    for name, steps in stages:
        print(f"\n=== {name} ===", flush=True)
        for script, args in steps:
            path = SCRIPTS / script
            if not path.is_file():
                print(f"ERROR: missing script {path}", file=sys.stderr)
                return 1
            print(f"+ python scripts/{script} {' '.join(args)}".rstrip(), flush=True)
            old_argv = sys.argv
            old_cwd = os.getcwd()
            sys.argv = [str(path), *args]
            try:
                os.chdir(ROOT)
                runpy.run_path(str(path), run_name="__main__")
            except SystemExit as exc:
                code = exc.code if isinstance(exc.code, int) else (0 if exc.code is None else 1)
                if code:
                    print(f"ERROR: {script} exited with {code}", file=sys.stderr)
                    return code
            finally:
                sys.argv = old_argv
                os.chdir(old_cwd)
    print("\nPipeline finished.")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true", help="list stages and scripts")
    ap.add_argument("--dry-run", action="store_true", help="print commands without running")
    ap.add_argument("--stage", help="run a single stage (name or number)")
    ap.add_argument("--from", dest="start", help="first stage (name or number)")
    ap.add_argument("--to", dest="stop", help="last stage (name or number)")
    args = ap.parse_args()

    if args.list:
        for i, (name, steps) in enumerate(STAGES, 1):
            print(f"{i}. {name}  ({len(steps)} scripts)")
        print("\nAlternative acquisition path (not run by default): "
              "scripts/02b_download_ena.py.")
        return 0

    stages = _select(args.start, args.stop, args.stage)
    if args.dry_run:
        _print_plan(stages)
        return 0
    return _run(stages)


if __name__ == "__main__":
    raise SystemExit(main())
