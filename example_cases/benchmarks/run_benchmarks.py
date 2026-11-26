#!/usr/bin/env python3
"""Runs the benchmarks in this folder with a solver companion and judges their checks, the
same way the `real_benchmarks_meet_their_references` test does, without building SolidCraft.

    python3 example_cases/benchmarks/run_benchmarks.py /path/to/cfd-solve OUTPUT_DIR [NAME...]

Writes OUTPUT_DIR/benchmarks.json (what each check measured) and prints a table.
"""
import csv
import json
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent


def column(path, name):
    with open(path, newline="") as f:
        rows = list(csv.DictReader(f))
    return [float(r[name]) for r in rows]


def value_of(spec, out, metrics):
    if "metric" in spec:
        v = metrics.get(spec["metric"])
        if v is None:
            return None
    elif "json" in spec:
        path = out / spec["json"]
        if not path.exists():
            return None
        v = json.loads(path.read_text())
        for key in spec["path"]:
            v = v[key]
    else:
        path = out / spec["file"]
        if not path.exists():
            return None
        ys = column(path, spec["column"])
        if "row" in spec:
            v = {"first": ys[0], "last": ys[-1], "max": max(ys), "min": min(ys)}[spec["row"]]
        else:
            xs = column(path, spec["at"]["column"])
            x = spec["at"]["value"]
            pairs = sorted(zip(xs, ys))
            v = None
            for (x0, y0), (x1, y1) in zip(pairs, pairs[1:]):
                if x0 <= x <= x1:
                    v = y0 if x1 == x0 else y0 + (y1 - y0) * (x - x0) / (x1 - x0)
                    break
            if v is None:
                return None
    if "over" in spec:
        under = value_of(spec["over"], out, metrics)
        if under is None:
            return None
        v /= under
    if "add" in spec:
        v += spec["add"]
    if "divide_by" in spec:
        v /= spec["divide_by"]
    return v


def judge(check, v):
    if "expect" in check:
        e = check["expect"]
        if not abs(v - e) <= check.get("tolerance", 0.05) * abs(e):
            return False
    if "at_least" in check and not v >= check["at_least"]:
        return False
    if "at_most" in check and not v <= check["at_most"]:
        return False
    return True


def main():
    solver, root = Path(sys.argv[1]), Path(sys.argv[2])
    only = sys.argv[3:]
    root.mkdir(parents=True, exist_ok=True)
    results = {}
    for path in sorted(HERE.glob("*.json")):
        name = path.stem
        if only and not any(name.startswith(o) for o in only):
            continue
        bench = json.loads(path.read_text())
        out = root / name
        if out.exists():
            for f in sorted(out.rglob("*"), reverse=True):
                f.unlink() if f.is_file() else f.rmdir()
            out.rmdir()
        case = root / f"{name}.case.json"
        case.write_text(json.dumps(bench["case"]))
        started = time.time()
        run = subprocess.run([str(solver), "run", str(case), "--out", str(out)], capture_output=True, text=True)
        seconds = time.time() - started
        if run.returncode != 0:
            results[name] = {"title": bench["title"], "error": (run.stderr or run.stdout)[-400:]}
            print(f"{name}: FAILED TO RUN: {(run.stderr or run.stdout)[-300:]}")
            continue
        metrics = json.loads((out / "manifest.json").read_text()).get("metrics", {})
        checks = []
        for check in bench["checks"]:
            v = value_of(check["value"], out, metrics)
            ok = v is not None and judge(check, v)
            off = (v - check["expect"]) / check["expect"] if v is not None and "expect" in check and check["expect"] else None
            checks.append({"what": check["what"], "value": v, "expect": check.get("expect"), "off": off, "pass": ok})
            print(f"{name}: {check['what']}: {v!r} vs {check.get('expect')!r} ({'' if off is None else f'{off*100:+.2f} %'}) {'pass' if ok else 'FAIL'}")
        results[name] = {"title": bench["title"], "domain": bench["domain"], "seconds": seconds, "checks": checks}
    (root / "benchmarks.json").write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n")


if __name__ == "__main__":
    main()
