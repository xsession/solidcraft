#!/usr/bin/env python3
"""Writes results.md from a run's benchmarks.json (from `real_benchmarks_meet_their_references`
or run_benchmarks.py):

    python3 example_cases/benchmarks/make_results.py OUTPUT_DIR/benchmarks.json "cfd-solvers-rs 0.18.0, 27 November 2026"
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent


def number(v):
    return "–" if v is None else f"{v:.4g}"


def main():
    runs = json.loads(Path(sys.argv[1]).read_text())
    label = sys.argv[2] if len(sys.argv) > 2 else ""
    lines = [
        "# Results of the last run",
        "",
        f"Every benchmark in this folder, run through `real_benchmarks_meet_their_references` ({label}).",
        "Each check's value, what it was held against, and whether it passed.",
        "",
        "| Benchmark | Domain | Check | Value | Expected | Off by | Pass | Run (s) |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for name in sorted(runs):
        run = runs[name]
        for i, check in enumerate(run.get("checks", [])):
            value, expect = check.get("value"), check.get("expect")
            off = f"{(value - expect) / abs(expect) * 100:+.2f} %" if value is not None and expect else "bound"
            first = i == 0
            lines.append(
                f"| {run['title'] if first else ''} | {run.get('domain', '') if first else ''} | {check['what']} | {number(value)} | "
                f"{number(expect) if expect is not None else 'bound'} | {off} | {'yes' if check['pass'] else '**no**'} | "
                f"{round(run['seconds']) if first else ''} |"
            )
    (HERE / "results.md").write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
