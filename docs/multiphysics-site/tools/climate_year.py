"""Run a Multiphysics study's climate year: one study per condition, every one run, the year
added up and the dust drawn where it settles.

The design holds a study of the model with a climate (multiphysics.model.climate) and,
usually, dust (multiphysics.model.dust). This script drives SolidCraft's command-line tool and
the solver companion; it reads nothing but their outputs:

    python3 climate_year.py prepare design.sc "Study name" work/
        makes one study per condition (multiphysics.model.climate_studies), saves the design as
        work/design.sc and writes each condition's case to work/cases/<n>.json
    python3 climate_year.py run work/ --jobs 4 [--dt 0.02 --end-time 30]
        runs every case with cfd-solve into work/runs/<n>/ (skips those already run). In the
        app a run is kept short (at most 50 million cell steps, a cautious time step); outside
        it --dt and --end-time set the time step and the air's time instead, and the run still
        stops once its goals settle. Check the result against a finer step before trusting it.
    python3 climate_year.py summary work/
        adds the year up (multiphysics.result.climate) into work/summary/ and draws
        work/summary/deposits-top.png and deposits-side.png

Set SOLIDCRAFT_CLI and SOLIDCRAFT_CFD_SOLVE to the two programs (default: on the PATH).
"""
import json
import math
import os
import struct
import subprocess
import sys
import zlib
from concurrent.futures import ThreadPoolExecutor

CLI = os.environ.get("SOLIDCRAFT_CLI", "solidcraft-cli")
SOLVER = os.environ.get("SOLIDCRAFT_CFD_SOLVE", "cfd-solve")


def cli(design, commands, save=None):
    """Run commands on a design with SolidCraft's command-line tool; their results."""
    script = os.path.abspath(os.path.join(os.path.dirname(design), ".climate-script.json"))
    out = script + ".out"
    with open(script, "w") as f:
        json.dump({"commands": [{"command": c, "params": p} for c, p in commands]}, f)
    args = [CLI, "run", script, "--in", design, "--results", "--json-out", out, "--quiet"]
    if save:
        args += ["--save", save]
    done = subprocess.run(args, capture_output=True, text=True)
    if done.returncode != 0:
        raise SystemExit(f"{CLI} failed: {done.stderr.strip()[-2000:]}")
    with open(out) as f:
        results = json.load(f)
    os.remove(script)
    os.remove(out)
    return results


def studies(design):
    return cli(design, [("multiphysics.study.list", {})])[0]["project"]["studies"]


def prepare(design, name, work):
    os.makedirs(os.path.join(work, "cases"), exist_ok=True)
    source = next((s for s in studies(design) if s["name"] == name), None)
    if source is None:
        raise SystemExit(f"the design has no study called {name}")
    saved = os.path.abspath(os.path.join(work, "design.sc"))
    made = cli(design, [("multiphysics.model.climate_studies", {"id": source["id"]})], save=saved)[0]
    exports = []
    for n, (study, condition) in enumerate(zip(made["studies"], made["conditions"]), 1):
        path = os.path.abspath(os.path.join(work, "cases", f"{n:02d}.json"))
        if os.path.exists(path):
            os.remove(path)
        exports.append(("multiphysics.case.export", {"id": study, "path": path}))
    cli(saved, exports)
    with open(os.path.join(work, "conditions.json"), "w") as f:
        json.dump({"source": source["id"], "studies": made["studies"], "conditions": made["conditions"]}, f, indent=1)
    print(f"{len(exports)} condition studies; cases in {os.path.join(work, 'cases')}")


def run(work, jobs, dt=None, end_time=None):
    with open(os.path.join(work, "conditions.json")) as f:
        plan = json.load(f)
    os.makedirs(os.path.join(work, "runs"), exist_ok=True)

    def one(n):
        case = os.path.join(work, "cases", f"{n:02d}.json")
        out = os.path.join(work, "runs", f"{n:02d}")
        if dt or end_time:
            with open(case) as f:
                c = json.load(f)
            step = dt or c["config"]["dt"]
            steps = max(1, math.ceil((end_time or step * c["output"]["steps"]) / step))
            c["config"]["dt"] = step
            c["output"] = {"steps": steps, "output_interval": steps}
            if "goal_control" in c:
                window = max(10, min(500, steps // 8))
                c["goal_control"].update({"window_steps": window, "min_steps": max(2 * window, steps // 4)})
            case = os.path.join(work, "cases", f"{n:02d}-run.json")
            with open(case, "w") as f:
                json.dump(c, f)
        if os.path.exists(os.path.join(out, "manifest.json")):
            return f"{n:02d} already run"
        if os.path.isdir(out):
            for name in os.listdir(out):
                os.remove(os.path.join(out, name))
            os.rmdir(out)
        with open(out + ".log", "w") as log:
            code = subprocess.run([SOLVER, "run", case, "--out", out], stdout=log, stderr=subprocess.STDOUT).returncode
        return f"{n:02d} {plan['conditions'][n - 1]}: {'done' if code == 0 else 'FAILED, see ' + out + '.log'}"

    with ThreadPoolExecutor(max_workers=jobs) as pool:
        for line in pool.map(one, range(1, len(plan["studies"]) + 1)):
            print(line, flush=True)


def summary(work):
    with open(os.path.join(work, "conditions.json")) as f:
        plan = json.load(f)
    folder = os.path.abspath(os.path.join(work, "summary"))
    if os.path.isdir(folder):
        for name in os.listdir(folder):
            os.remove(os.path.join(folder, name))
        os.rmdir(folder)
    runs = [
        {"study": study, "directory": os.path.abspath(os.path.join(work, "runs", f"{n:02d}"))}
        for n, study in enumerate(plan["studies"], 1)
        if os.path.exists(os.path.join(work, "runs", f"{n:02d}", "manifest.json"))
    ]
    result = cli(os.path.join(work, "design.sc"), [("multiphysics.result.climate", {"id": plan["source"], "runs": runs, "directory": folder})])[0]
    print(f"dust in: {result['dust_in_g_per_year']:.4f} g a year, settled: {result['deposited_g_per_year']:.4f} g")
    for p in result["parts"][:12]:
        print(f"  {p['body']:<14} {p['dust_g_per_year']:.5f} g/yr  humid dust {p['humid_dust_hours']:5.0f} h  below dew point {p['condensing_hours']:5.0f} h")
    draw_maps(folder)


# --- Pictures, with nothing but the standard library -----------------------------------------
def png(path, pixels, width, height):
    """Write RGB rows (bytes of width*3 each) as a PNG."""
    raw = b"".join(b"\x00" + row for row in pixels)

    def chunk(kind, data):
        return struct.pack(">I", len(data)) + kind + data + struct.pack(">I", zlib.crc32(kind + data) & 0xFFFFFFFF)

    with open(path, "wb") as f:
        f.write(b"\x89PNG\r\n\x1a\n")
        f.write(chunk(b"IHDR", struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)))
        f.write(chunk(b"IDAT", zlib.compress(raw, 9)))
        f.write(chunk(b"IEND", b""))


def colour(t):
    """A light-to-dark sequential ramp: pale yellow through orange to dark brown."""
    stops = [(0.0, (255, 247, 188)), (0.35, (254, 196, 79)), (0.7, (217, 95, 14)), (1.0, (102, 37, 6))]
    t = min(max(t, 0.0), 1.0)
    for (a, ca), (b, cb) in zip(stops, stops[1:]):
        if t <= b:
            f = (t - a) / (b - a)
            return tuple(int(ca[k] + f * (cb[k] - ca[k])) for k in range(3))
    return stops[-1][1]


def draw_maps(folder):
    rows = []
    with open(os.path.join(folder, "climate-deposits.csv")) as f:
        next(f)
        for line in f:
            x, y, z, part, d, g, humid, _ = line.rstrip("\n").split(",", 7)
            rows.append((float(x), float(y), float(z), part, float(g), float(humid)))
    if not rows:
        print("no deposits to draw")
        return
    for view, (a, b) in {"top": (0, 1), "side": (0, 2)}.items():
        lo = [min(r[a] for r in rows), min(r[b] for r in rows)]
        hi = [max(r[a] for r in rows), max(r[b] for r in rows)]
        mm = 1.0  # a pixel per bin of 1 mm, drawn 4 px wide
        nx, ny = int((hi[0] - lo[0]) / mm) + 1, int((hi[1] - lo[1]) / mm) + 1
        grid = [[0.0] * nx for _ in range(ny)]
        for r in rows:
            i, j = int((r[a] - lo[0]) / mm), int((r[b] - lo[1]) / mm)
            grid[j][i] += r[4] / (mm * mm) * 1e6  # g/m² a year
        peak = max(max(row) for row in grid) or 1.0
        scale = 4
        pixels = []
        for j in reversed(range(ny)):
            line = bytearray()
            for i in range(nx):
                v = grid[j][i]
                rgb = (250, 250, 250) if v <= 0 else colour(math.log10(1 + 99 * v / peak) / 2)
                line += bytes(rgb) * scale
            for _ in range(scale):
                pixels.append(bytes(line))
        path = os.path.join(folder, f"deposits-{view}.png")
        png(path, pixels, nx * scale, ny * scale)
        print(f"{path}: {nx} x {ny} mm, darkest {peak:.3g} g/m² a year")


def main(argv):
    if len(argv) >= 4 and argv[0] == "prepare":
        prepare(argv[1], argv[2], argv[3])
    elif len(argv) >= 2 and argv[0] == "run":
        option = lambda name, kind: kind(argv[argv.index(name) + 1]) if name in argv else None
        run(argv[1], option("--jobs", int) or 2, option("--dt", float), option("--end-time", float))
    elif len(argv) >= 2 and argv[0] == "summary":
        summary(argv[1])
    else:
        raise SystemExit(__doc__)


if __name__ == "__main__":
    main(sys.argv[1:])
