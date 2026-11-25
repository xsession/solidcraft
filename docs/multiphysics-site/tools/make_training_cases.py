"""Build the training cases in SolidCraft, run them for real and write their results.

Each case is drawn here from scratch (our own geometry and settings), set up as a study of the
model, run with the native companion, and checked against an exact or textbook reference where
one exists. The script drives a running SolidCraft through its control channel (commands only,
no clicks), like make_study.py. It writes:

- one design per case into OUTDIR (kept out of the repository: they hold the runs);
- OUTDIR/results.json: every case's numbers, which the Training cases page quotes;
- OUTDIR/<case>.png: a picture of each case's result.

    SOLIDCRAFT_CFD_SOLVE=/path/to/cfd-solve solidcraft --control 7979 &
    SOLIDCRAFT_CONTROL_TOKEN=$(cat ~/.local/share/solidcraft/control/7979.token) \\
        python3 make_training_cases.py OUTDIR [case …]
"""
import functools
import json
import math
import os
import socket
import sys
import time
import zipfile

PORT = int(os.environ.get("PORT", "7979"))
TOKEN = os.environ["SOLIDCRAFT_CONTROL_TOKEN"]
OUTDIR = os.path.abspath(sys.argv[1])
ONLY = set(sys.argv[2:])
os.makedirs(OUTDIR, exist_ok=True)

print = functools.partial(print, flush=True)  # noqa: A001 (progress shows as it happens)
sock = socket.create_connection(("127.0.0.1", PORT), timeout=3600)
reader = sock.makefile("r")
counter = 0

AIR_NU = 1.5e-5  # m²/s, air at 20 °C
AIR_RHO = 1.2  # kg/m³
AIR_CP = 1005.0  # J/(kg·K)


def call(method, params=None):
    global counter
    counter += 1
    sock.sendall((json.dumps({"id": counter, "method": method, "params": params or {}, "token": TOKEN}) + "\n").encode())
    reply = json.loads(reader.readline())
    # A command that outlasts the channel's wait still runs in the app: ask for its result.
    while not reply.get("ok", False) and reply.get("running"):
        time.sleep(5)
        counter += 1
        sock.sendall((json.dumps({"id": counter, "method": "control.pending", "params": {}, "token": TOKEN}) + "\n").encode())
        pending = json.loads(reader.readline())
        result = pending.get("result") or {}
        if result.get("done"):
            reply = result.get("reply") or {}
        elif result.get("pending") is None:
            raise RuntimeError(f"{method}: the app lost the request")
    if not reply.get("ok", False):
        raise RuntimeError(f"{method} failed: {json.dumps(reply)[:900]}")
    return reply["result"]


def run(command, params=None):
    return call("engine.execute", {"command": command, "params": params or {}})


def new_design():
    run("file.new", {})


class Study:
    """A study of the model, built from the design's bodies."""

    def __init__(self, name, cell_mm, flow, analysis="internal", cut_cells=False):
        params = {"name": name, "analysis": analysis, "flow": flow}
        if cell_mm:
            params["cell_mm"] = cell_mm
        reply = run("multiphysics.model.create", params)
        self.id = reply["study"]["id"]
        self.view = reply["study"]
        self.notes = list(reply.get("notes", []))
        if cut_cells:
            self.set(cut_cells=True)

    def revision(self):
        return next(s for s in run("multiphysics.study.list")["project"]["studies"] if s["id"] == self.id)["revision"]

    def do(self, command, **params):
        reply = run(command, {"id": self.id, "expected_revision": self.revision(), **params})
        self.notes = list(reply.get("notes", self.notes))
        if "study" in reply:
            self.view = reply["study"]
        return reply

    def set(self, **params):
        return self.do("multiphysics.model.set", **params)

    def condition(self, name, body, point, role):
        return self.do("multiphysics.model.condition", name=name, face={"body": body, "point": point}, role=role)

    def goal(self, name, kind, **extra):
        return self.do("multiphysics.model.goal", name=name, kind=kind, **extra)

    def dimensions(self):
        return self.view["case"]["mesh"]["dimensions"]

    def solve(self, label):
        """Run, wait, keep; returns the kept run's id and its manifest metrics."""
        started = time.time()
        job = run("multiphysics.run.start", {"id": self.id})
        job_id = job.get("job_id") or job.get("job", {}).get("id") or job.get("id")
        while True:
            time.sleep(2)
            status = run("multiphysics.run.status")
            jobs = status.get("jobs", status if isinstance(status, list) else [])
            mine = next((j for j in jobs if j.get("id") == job_id), None) or {}
            state = mine.get("state")
            if state in ("completed", "succeeded", "done", "finished") or mine.get("evidence"):
                break
            if state in ("failed", "cancelled", "error"):
                raise RuntimeError(f"run {label} {state}: {json.dumps(mine)[:900]}")
            if time.time() - started > 3 * 3600:
                raise RuntimeError(f"run {label} took over three hours")
        kept = run("multiphysics.run.keep", {"job_id": job_id})
        run_id = kept["study"]["runs"][-1]["id"]
        seconds = time.time() - started
        print(f"  {label}: run {run_id} in {seconds:.0f} s")
        return run_id, seconds

    def query(self, run_id, query):
        return run("multiphysics.result.query", {"id": self.id, "run_id": run_id, "query": query})



def saved_run(path, study_id, run_id):
    """The kept run as the saved design holds it: its manifest (metrics) and artifact names."""
    with zipfile.ZipFile(path) as z:
        document = json.loads(z.read("document.json"))
    project = document.get("multiphysics") or document.get("document", {}).get("multiphysics")
    study = next(s for s in project["studies"] if s["id"] == study_id)
    record = next(r for r in study["runs"] if r["id"] == run_id)
    return {"metrics": record["manifest"].get("metrics", {}), "artifacts": [a["name"] for a in record["artifacts"]], "case": record["input_case"]}


def last_flow(record):
    return sorted(n for n in record["artifacts"] if n.startswith("flow_") and n.endswith(".vtu"))[-1]


def probe(study, run_id, artifact, field, point_m):
    reply = study.query(run_id, {"kind": "probe", "artifact": artifact, "field": field, "point": point_m})
    return reply["result"]


def picture(name):
    path = os.path.join(OUTDIR, f"{name}.png")
    try:
        call("ui.screenshot", {"path": path})
    except RuntimeError as error:
        print("  no picture:", error)
        return None
    return path


# Earlier runs' results stay: a run of some cases adds to them.
try:
    with open(os.path.join(OUTDIR, "results.json")) as previous:
        RESULTS = json.load(previous)
except (OSError, ValueError):
    RESULTS = {}


def case(fn):
    key = fn.__name__
    if ONLY and key not in ONLY:
        return fn
    print(f"== {key}")
    started = time.time()
    try:
        RESULTS[key] = fn()
    except Exception as error:  # keep going: one failed case is a result too
        RESULTS[key] = {"error": str(error)}
        print("  failed:", error)
    RESULTS[key]["wall_seconds"] = round(time.time() - started)
    with open(os.path.join(OUTDIR, "results.json"), "w") as out:
        json.dump(RESULTS, out, indent=2)
    return fn


def icosphere(name, center, diameter, level=4):
    """A sphere as a faceted solid: an icosahedron split `level` times (20·4^level triangles,
    within 0.3 % of the radius at level 4), written as STL, inserted and made a solid. The
    kernel's own tessellation of a sphere is far finer than a study's cells need and passes the
    case's size limit."""
    t = (1 + 5 ** 0.5) / 2
    verts = [(-1, t, 0), (1, t, 0), (-1, -t, 0), (1, -t, 0), (0, -1, t), (0, 1, t), (0, -1, -t), (0, 1, -t), (t, 0, -1), (t, 0, 1), (-t, 0, -1), (-t, 0, 1)]
    verts = [tuple(c / math.sqrt(sum(x * x for x in v)) for c in v) for v in verts]
    faces = [(0, 11, 5), (0, 5, 1), (0, 1, 7), (0, 7, 10), (0, 10, 11), (1, 5, 9), (5, 11, 4), (11, 10, 2), (10, 7, 6), (7, 1, 8),
             (3, 9, 4), (3, 4, 2), (3, 2, 6), (3, 6, 8), (3, 8, 9), (4, 9, 5), (2, 4, 11), (6, 2, 10), (8, 6, 7), (9, 8, 1)]
    for _ in range(level):
        cache, out = {}, []

        def mid(a, b):
            key = (min(a, b), max(a, b))
            if key not in cache:
                m = [(verts[a][k] + verts[b][k]) / 2 for k in range(3)]
                r = math.sqrt(sum(x * x for x in m))
                verts.append(tuple(x / r for x in m))
                cache[key] = len(verts) - 1
            return cache[key]

        for a, b, c in faces:
            ab, bc, ca = mid(a, b), mid(b, c), mid(c, a)
            out += [(a, ab, ca), (b, bc, ab), (c, ca, bc), (ab, bc, ca)]
        faces = out
    r = diameter / 2
    path = os.path.join(OUTDIR, f"{name}.stl")
    with open(path, "w") as f:
        f.write(f"solid {name}\n")
        for tri in faces:
            f.write("facet normal 0 0 0\nouter loop\n")
            for i in tri:
                f.write("vertex " + " ".join(f"{center[k] + r * verts[i][k]:.6f}" for k in range(3)) + "\n")
            f.write("endloop\nendfacet\n")
        f.write(f"endsolid {name}\n")
    return solid_from_stl(path, name)


def closed_stl(paths, out_path):
    """Join STL surfaces (a channel and its opening caps) into one closed, consistently oriented
    surface: each cap is flipped where that makes every edge appear once in each direction."""
    import itertools

    import numpy as np

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from stl_parts import load

    parts = [load(p) for p in paths]
    span = max(np.ptp(np.concatenate(parts).reshape(-1, 3), axis=0).max(), 1.0)

    def edges(tris):
        keys = np.round(tris / (span * 1e-7)).astype(np.int64)
        out = {}
        for t in keys:
            for a, b in ((0, 1), (1, 2), (2, 0)):
                e = (tuple(t[a]), tuple(t[b]))
                out[e] = out.get(e, 0) + 1
        return out

    best = None
    for flips in itertools.product([False, True], repeat=len(parts) - 1):
        tris = np.concatenate([parts[0]] + [p[:, ::-1] if f else p for p, f in zip(parts[1:], flips)])
        e = edges(tris)
        bad = sum(1 for (a, b), n in e.items() if n != 1 or e.get((b, a), 0) != 1)
        if best is None or bad < best[0]:
            best = (bad, tris)
    with open(out_path, "w") as f:
        f.write("solid joined\n")
        for t in best[1]:
            f.write("facet normal 0 0 0\nouter loop\n" + "".join(f"vertex {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n" for v in t) + "endloop\nendfacet\n")
        f.write("endsolid joined\n")
    return best[0]


def solid_from_stl(path, name):
    """Insert an STL and make it a solid; returns the body's name."""
    body = run("file.insert_mesh", {"path": path, "units": "mm", "name": name})["created_bodies"][0]
    run("mesh.convert", {"body": body, "name": f"{name} solid"})
    return body


def key_of(name):
    return name.lower().replace(", ", "_").replace(" ", "_")


def save(name):
    path = os.path.join(OUTDIR, f"{name}.sc")
    run("file.save_as", {"path": path})
    return path


# ---------------------------------------------------------------------------------------------
# Laminar flow in straight passages: the pressure gradient where the flow is fully developed,
# against the exact laminar value  dp/dx = (f·Re) · ν · U / (2 · Dh²)  (kinematic pressure),
# with f·Re = 64 for a round pipe, 56.91 for a square duct and 82.34 for a 1 : 8 rectangle
# (Shah & London, "Laminar Flow Forced Convection in Ducts", 1978).


def passage(name, bore, length, wall, cell, speed, f_re, dh_mm, measure, end_time, cut_cells=False, round_bore=None):
    new_design()
    outer = [bore[0] + 2 * wall, bore[1] + 2 * wall]
    run("solid.box", {"length": length, "width": outer[0], "height": outer[1], "corner": [0, 0, 0], "operation": "new", "body_name": "Wall"})
    if round_bore:
        run("solid.cylinder", {"diameter": round_bore, "height": length + 10, "base": [-5, outer[0] / 2, outer[1] / 2], "axis": "X", "operation": "cut"})
    else:
        run("solid.box", {"length": length + 10, "width": bore[0], "height": bore[1], "corner": [-5, wall, wall], "operation": "cut"})
    flow = {"kinematic_viscosity_m2_per_s": AIR_NU, "density_kg_per_m3": AIR_RHO, "end_time_s": end_time, "converge_tolerance": 2e-4}
    study = Study(name, cell, flow, cut_cells=cut_cells)
    # Inlet and outlet on the wall's end faces (a point on the end face's solid rim).
    study.condition("Inlet", "Wall", [0, wall / 2, outer[1] / 2], {"kind": "inlet", "speed_m_per_s": speed})
    study.condition("Outlet", "Wall", [length, wall / 2, outer[1] / 2], {"kind": "outlet", "pressure_m2_per_s2": 0.0})
    mid = [outer[0] / 2, outer[1] / 2]
    for i, x in enumerate(measure):
        study.goal(f"p{i}", "point_pressure", point_mm=[x, mid[0], mid[1]])
    study.goal("q_out", "face_volume_flow", face={"body": "Wall", "point": [length, wall / 2, outer[1] / 2]})
    run_id, seconds = study.solve(name)
    dims = study.dimensions()
    picture(key_of(name))
    m = saved_run(save(key_of(name)), study.id, run_id)["metrics"]
    # The flow through the bore over its true area: the mean speed the exact value is for.
    area_m2 = (math.pi * round_bore**2 / 4 if round_bore else bore[0] * bore[1]) * 1e-6
    mean = abs(m["goal_q_out"]) / area_m2
    expected = f_re * AIR_NU * mean / (2 * (dh_mm / 1000.0) ** 2)
    # Point pressure goals are in pascals (kinematic pressure × density).
    measured = (m["goal_p0"] - m["goal_p1"]) / ((measure[1] - measure[0]) / 1000.0) / AIR_RHO
    return {
        "cells": dims,
        "cell_count": dims[0] * dims[1] * dims[2],
        "cell_mm": cell,
        "cut_cells": cut_cells,
        "mean_speed_m_per_s": mean,
        "reynolds": mean * dh_mm / 1000.0 / AIR_NU,
        "dpdx_measured_m2_per_s2_per_m": measured,
        "dpdx_exact_m2_per_s2_per_m": expected,
        "error_percent": (measured / expected - 1) * 100,
        "goals_converged": m.get("goals_converged"),
        "steps": m.get("completed_steps") or m.get("steps"),
        "simulated_s": m.get("time_s"),
        "run_seconds": seconds,
        "continuity_l2": m.get("continuity_l2"),
        "notes": study.notes,
    }


@case
def round_pipe():
    # Bore 16 mm, 120 mm long, 1.5 mm cells; uniform inflow 0.05 m/s (Re ≈ 53, developed after
    # about 0.06·Re·D ≈ 50 mm); measured from 70 to 110 mm.
    return passage("Round pipe", [16, 16], 120, 4, 1.5, 0.05, 64.0, 16, [70, 110], 4.0, round_bore=16)


@case
def round_pipe_cut_cells():
    return passage("Round pipe, cut cells", [16, 16], 120, 4, 1.5, 0.05, 64.0, 16, [70, 110], 4.0, cut_cells=True, round_bore=16)


@case
def square_duct():
    # 8 × 8 mm, 60 mm long, 0.8 mm cells (10 across): the grid holds it exactly.
    return passage("Square duct", [8, 8], 60, 1.6, 0.8, 0.05, 56.91, 8, [30, 52], 3.0)


@case
def thin_gap():
    # A 2 × 16 mm slot (1 : 8), 30 mm long, 0.4 mm cells (5 across the gap).
    dh = 2 * 2 * 16 / (2 + 16)
    return passage("Thin gap", [16, 2], 30, 0.8, 0.4, 0.05, 82.34, dh, [16, 27], 0.85)


# ---------------------------------------------------------------------------------------------
@case
def sphere_in_channel():
    # A 10 mm sphere on the axis of a 30 × 30 mm channel, 90 mm long, 1.7 mm cells, 0.1 m/s
    # (Re = 67 on the sphere). No exact answer: the wake's length behind it and the extra
    # pressure drop it causes.
    new_design()
    run("solid.box", {"length": 90, "width": 34, "height": 34, "corner": [0, 0, 0], "operation": "new", "body_name": "Channel"})
    run("solid.box", {"length": 100, "width": 30, "height": 30, "corner": [-5, 2, 2], "operation": "cut"})
    icosphere("Sphere", [30, 17, 17], 10)
    flow = {"kinematic_viscosity_m2_per_s": AIR_NU, "density_kg_per_m3": AIR_RHO, "end_time_s": 1.5, "converge_tolerance": 5e-4}
    study = Study("Sphere in channel", 1.7, flow)
    study.condition("Inlet", "Channel", [0, 1, 17], {"kind": "inlet", "speed_m_per_s": 0.1})
    study.condition("Outlet", "Channel", [90, 1, 17], {"kind": "outlet", "pressure_m2_per_s2": 0.0})
    study.goal("p_in", "point_pressure", point_mm=[5, 17, 17])
    study.goal("max_speed", "global_max_speed")
    run_id, seconds = study.solve("sphere")
    picture("sphere_in_channel")
    record = saved_run(save("sphere_in_channel"), study.id, run_id)
    m = record["metrics"]
    # Along the axis behind the sphere (result coordinates: m from the domain's corner, which
    # is the design's origin here): where the flow along x turns forward again.
    artifact = last_flow(record)
    wake_end_mm = 35.0
    for x10 in range(352, 800, 5):
        u = probe(study, run_id, artifact, "velocity_m_per_s", [x10 / 10000.0, 0.017, 0.017])["value"]
        if u and u[0] < 0:
            wake_end_mm = x10 / 10.0
    return {
        "reynolds": 0.1 * 0.010 / AIR_NU,
        "wake_length_over_d": (wake_end_mm - 35.0) / 10.0,
        "inlet_pressure_pa": m.get("goal_p_in"),
        "max_speed_m_per_s": m.get("goal_max_speed"),
        "goals_converged": m.get("goals_converged"),
        "run_seconds": seconds,
        "notes": study.notes,
    }


# ---------------------------------------------------------------------------------------------
@case
def battery_pack():
    # Six 46 × 80 mm cylindrical cells (2 × 3) standing in an air duct, 5 W each; air at 0.2 m/s
    # (laminar: the solver has no turbulence model in 3-D flow). The check: the energy the air
    # carries out equals the 30 W the cells release, and the outlet's mean rise is
    # P / (ρ·cp·Q).
    new_design()
    run("solid.box", {"length": 300, "width": 150, "height": 100, "corner": [0, 0, 0], "operation": "new", "body_name": "Duct"})
    run("solid.box", {"length": 310, "width": 142, "height": 92, "corner": [-5, 4, 4], "operation": "cut"})
    names = []
    for i, x in enumerate([90, 150, 210]):
        for j, y in enumerate([45, 105]):
            name = f"Cell{i * 2 + j + 1}"
            names.append(name)
            run("solid.cylinder", {"diameter": 46, "height": 80, "base": [x, y, 8], "operation": "new", "body_name": name})
    flow = {"kinematic_viscosity_m2_per_s": AIR_NU, "density_kg_per_m3": AIR_RHO, "end_time_s": 4.0, "converge_tolerance": 1e-3}
    study = Study("Battery pack", 6.0, flow)
    study.condition("Inlet", "Duct", [0, 2, 50], {"kind": "inlet", "speed_m_per_s": 0.2})
    study.condition("Outlet", "Duct", [300, 2, 50], {"kind": "outlet", "pressure_m2_per_s2": 0.0})
    study.do("multiphysics.model.heat", enabled=True, inflow_temperature_k=298.15)
    study.do("multiphysics.model.material", bodies=names, conductivity_w_per_m_k=3.0, power_w=5.0, max_temperature_k=318.15)
    study.do("multiphysics.model.material", bodies=["Duct"], conductivity_w_per_m_k=0.2, power_w=0.0)
    study.goal("hottest", "global_max_solid_temperature")
    study.goal("outlet_flow", "face_volume_flow", face={"body": "Duct", "point": [300, 2, 50]})
    run_id, seconds = study.solve("battery")
    report = study.query(run_id, {"kind": "components"})["result"]["rows"]
    picture("battery_pack")
    m = saved_run(save("battery_pack"), study.id, run_id)["metrics"]
    q = abs(m.get("goal_outlet_flow") or 0.142 * 0.092 * 0.2)
    rise = 30.0 / (AIR_RHO * AIR_CP * q)
    return {
        "total_power_w": m.get("thermal_total_power_w"),
        "relative_energy_imbalance": m.get("thermal_relative_energy_imbalance"),
        "outlet_mean_rise_expected_k": rise,
        "cells": [{"name": r.get("name"), "max_c": r.get("max_temperature_k", 0) - 273.15, "status": r.get("status")} for r in report],
        "hottest_c": (m.get("goal_hottest") or 0) - 273.15,
        "reynolds_cell": 0.2 * 0.046 / AIR_NU,
        "goals_converged": m.get("goals_converged"),
        "run_seconds": seconds,
        "metrics": {k: v for k, v in m.items() if k.startswith("thermal_")},
        "notes": study.notes,
    }


# ---------------------------------------------------------------------------------------------
@case
def fan_in_duct():
    # A fan blowing into the 8 × 8 mm square duct, 60 mm long, against ambient pressure. Its
    # curve falls straight from 0.1 Pa at no flow to none at 1.2 L/min. The duct's laminar system
    # curve (developed flow) is Δp = ρ·56.91·ν·L·Q / (2·a²·A); where they cross is where it runs.
    new_design()
    wall, bore, length = 1.6, 8.0, 60.0
    outer = bore + 2 * wall
    run("solid.box", {"length": length, "width": outer, "height": outer, "corner": [0, 0, 0], "operation": "new", "body_name": "Wall"})
    run("solid.box", {"length": length + 10, "width": bore, "height": bore, "corner": [-5, wall, wall], "operation": "cut"})
    flow = {"kinematic_viscosity_m2_per_s": AIR_NU, "density_kg_per_m3": AIR_RHO, "end_time_s": 1.1, "converge_tolerance": 2e-4}
    study = Study("Fan in duct", 0.8, flow)
    p0, qmax = 0.1, 2e-5
    study.condition("Outlet", "Wall", [length, wall / 2, outer / 2], {"kind": "outlet", "pressure_m2_per_s2": 0.0})
    study.do("multiphysics.model.fan", face={"body": "Wall", "point": [0, wall / 2, outer / 2]}, curve=[[0.0, p0], [qmax, 0.0]], name="Fan")
    study.goal("flow", "face_volume_flow", face={"body": "Wall", "point": [length, wall / 2, outer / 2]})
    run_id, seconds = study.solve("fan")
    picture("fan_in_duct")
    m = saved_run(save("fan_in_duct"), study.id, run_id)["metrics"]
    a, area = bore / 1000.0, (bore / 1000.0) ** 2
    k = AIR_RHO * 56.91 * AIR_NU * (length / 1000.0) / (2 * a * a * area)
    q_expected = p0 / (k + p0 / qmax)
    q = abs(m.get("goal_flow") or m.get("boundary_outflow_m3_per_s") or 0.0)
    return {
        "flow_m3_per_s": q,
        "flow_expected_m3_per_s": q_expected,
        "error_percent": (q / q_expected - 1) * 100 if q_expected else None,
        "goals_converged": m.get("goals_converged"),
        "run_seconds": seconds,
        "notes": study.notes,
    }


# ---------------------------------------------------------------------------------------------
# The downloaded cases' own geometry, measured from their STL files (stl_parts.py) and drawn
# again here with exact shapes; nothing else is taken from those cases. GEOMETRY points at the
# extracted files (outside the repository).
GEOMETRY = os.environ.get("TRAINING_GEOMETRY", os.path.join(os.path.dirname(OUTDIR), "geometry"))


@case
def battery_cooling_exact():
    # 20 cylindrical cells, 46 mm × 82 mm, on a staggered 25 mm pitch (measured), in the case's
    # own air channel (its STL, made a solid and cut out of a block); air in through an 80 × 40 mm
    # opening at 0.5 m/s (laminar), 2 W per cell (the case has 100 W per cell, which needs
    # turbulent air at about 10 m/s: out of a laminar solver's reach). Check: the energy balance.
    new_design()
    stl = os.path.join(GEOMETRY, "batteryCooling/cases-12/batteryCooling/cad/stl")
    joined = os.path.join(OUTDIR, "battery_channel.stl")
    open_edges = closed_stl([os.path.join(stl, n) for n in ("channel.stl", "inlet.stl", "outlet.stl")], joined)
    print("  channel joined, edges left open:", open_edges)
    channel = solid_from_stl(joined, "Channel")
    run("solid.box", {"length": 325, "width": 211, "height": 106, "corner": [-50, -172.5, -12], "operation": "new", "body_name": "Wall"})
    run("solid.combine", {"target": "Wall", "tools": [channel], "operation": "cut"})
    names = []
    for i in range(10):
        x = i * 25.0
        rows = [-89.0, 0.0] if i % 2 == 0 else [-135.0, -44.0]
        for y in rows:
            name = f"Cell{len(names) + 1}"
            names.append(name)
            run("solid.cylinder", {"diameter": 46, "height": 82, "base": [x, y, 0], "operation": "new", "body_name": name})
    flow = {"kinematic_viscosity_m2_per_s": AIR_NU, "density_kg_per_m3": AIR_RHO, "end_time_s": 3.7, "converge_tolerance": 1e-3}
    study = Study("Battery cooling (case geometry)", 8.0, flow)
    # A point on the solid part of each end face (the openings are 80 × 40 mm at y −106…−26, z 20…60).
    study.condition("Inlet", "Wall", [-50, -160, 80], {"kind": "inlet", "speed_m_per_s": 0.5})
    study.condition("Outlet", "Wall", [275, -160, 80], {"kind": "outlet", "pressure_m2_per_s2": 0.0})
    study.do("multiphysics.model.heat", enabled=True, inflow_temperature_k=298.15)
    study.do("multiphysics.model.material", bodies=names, conductivity_w_per_m_k=3.0, power_w=2.0, max_temperature_k=333.15)
    study.do("multiphysics.model.material", bodies=["Wall"], conductivity_w_per_m_k=0.2, power_w=0.0)
    study.goal("hottest", "global_max_solid_temperature")
    study.goal("q_out", "face_volume_flow", face={"body": "Wall", "point": [275, -160, 80]})
    run_id, seconds = study.solve("battery (case geometry)")
    rows = study.query(run_id, {"kind": "components"})["result"]["rows"]
    picture("battery_cooling_exact")
    m = saved_run(save("battery_cooling_exact"), study.id, run_id)["metrics"]
    q = abs(m.get("goal_q_out") or 0.0)
    return {
        "cells": study.dimensions(),
        "inflow_m3_per_s": q,
        "air_rise_expected_k": 40.0 / (AIR_RHO * AIR_CP * q) if q else None,
        "total_power_w": m.get("thermal_total_power_w"),
        "relative_energy_imbalance": m.get("thermal_relative_energy_imbalance"),
        "hottest_c": (m.get("goal_hottest") or 0) - 273.15,
        "parts": [{"name": r["name"], "max_c": r["max_temperature_k"] - 273.15, "status": r["status"]} for r in rows],
        "goals_converged": m.get("goals_converged"),
        "simulated_s": m.get("time_s"),
        "run_seconds": seconds,
        "notes": study.notes,
    }


@case
def sphere_exact():
    # The case's square channel (100 × 100 mm, 500 mm long, along z) and 20 mm sphere 100 mm
    # from the inlet, measured; air at 0.075 m/s (Re = 100 on the sphere). Cells split twice as
    # fine around the sphere. Reference: the closed wake behind a sphere at Re = 100 reaches
    # about 0.9 diameters past its rear (Johnson & Patel, J. Fluid Mech. 378, 1999).
    new_design()
    run("solid.box", {"length": 104, "width": 104, "height": 500, "corner": [-52, -52, -100], "operation": "new", "body_name": "Channel"})
    run("solid.box", {"length": 100, "width": 100, "height": 510, "corner": [-50, -50, -105], "operation": "cut"})
    icosphere("Sphere", [0, 0, 0], 20)
    flow = {"kinematic_viscosity_m2_per_s": AIR_NU, "density_kg_per_m3": AIR_RHO, "end_time_s": 4.0, "converge_tolerance": 5e-4}
    study = Study("Sphere (case geometry)", 8.0, flow)
    study.condition("Inlet", "Channel", [-51, 0, -100], {"kind": "inlet", "speed_m_per_s": 0.075})
    study.condition("Outlet", "Channel", [-51, 0, 400], {"kind": "outlet", "pressure_m2_per_s2": 0.0})
    study.do("multiphysics.model.refine", name="around sphere", body="Sphere", factor=2, margin_cells=0)
    study.goal("behind", "point_pressure", point_mm=[0, 0, 30])
    run_id, seconds = study.solve("sphere (case geometry)")
    picture("sphere_exact")
    record = saved_run(save("sphere_exact"), study.id, run_id)
    m = record["metrics"]
    artifact = last_flow(record)
    # Result coordinates: m from the domain's low corner, (-52, -52, -100) mm here.
    to_m = lambda p: [(p[0] + 52) / 1000.0, (p[1] + 52) / 1000.0, (p[2] + 100) / 1000.0]
    wake_end = 10.0
    for z10 in range(102, 600, 5):
        u = probe(study, run_id, artifact, "velocity_m_per_s", to_m([0, 0, z10 / 10.0]))["value"]
        if u and u[2] < 0:
            wake_end = z10 / 10.0
    return {
        "cells": study.dimensions(),
        "reynolds": 0.075 * 0.020 / AIR_NU,
        "wake_length_over_d": (wake_end - 10.0) / 20.0,
        "wake_reference_over_d": 0.88,
        "goals_converged": m.get("goals_converged"),
        "simulated_s": m.get("time_s"),
        "run_seconds": seconds,
        "notes": study.notes,
    }


@case
def teg_stack_exact():
    # The TEG test device's parts, measured: a 100 × 150 × 20 mm cooler and heater, two
    # 25 × 25 × 35 mm copper blocks and the 25 × 25 × 2 mm thermoelectric layer, stacked
    # (the parts as listed; the case's own gap and insulation are not taken). 5 W in the heater,
    # the cooler's underside held at 20 °C, still air around. Check: the drop across the layer is
    # P·L/(k·A) = 5 × 0.002 / (1.5 × 0.025²) = 10.7 K when all the heat passes through it
    # (the air beside it carries about 1 %).
    new_design()
    run("solid.box", {"length": 100, "width": 150, "height": 20, "corner": [0, 0, 0], "operation": "new", "body_name": "Cooler"})
    run("solid.box", {"length": 25, "width": 25, "height": 35, "corner": [37.5, 62.5, 20], "operation": "new", "body_name": "CopperBlock"})
    run("solid.box", {"length": 25, "width": 25, "height": 2, "corner": [37.5, 62.5, 55], "operation": "new", "body_name": "TEG"})
    run("solid.box", {"length": 25, "width": 25, "height": 35, "corner": [37.5, 62.5, 57], "operation": "new", "body_name": "CopperBlock2"})
    run("solid.box", {"length": 100, "width": 150, "height": 20, "corner": [0, 0, 92], "operation": "new", "body_name": "Heater"})
    flow = {"kinematic_viscosity_m2_per_s": AIR_NU, "density_kg_per_m3": AIR_RHO, "end_time_s": 0.05, "converge_tolerance": None}
    study = Study("TEG stack (case parts)", 5.0, flow)
    study.do("multiphysics.model.heat", enabled=True, inflow_temperature_k=293.15)
    for bodies, k, power in [(["Cooler", "Heater"], 200.0, 0.0), (["CopperBlock", "CopperBlock2"], 390.0, 0.0), (["TEG"], 1.5, 0.0)]:
        study.do("multiphysics.model.material", bodies=bodies, conductivity_w_per_m_k=k, power_w=power)
    study.do("multiphysics.model.material", bodies=["Heater"], conductivity_w_per_m_k=200.0, power_w=5.0)
    study.do("multiphysics.model.thermal_face", name="Cold side", face={"body": "Cooler", "point": [50, 75, 0]},
             condition={"kind": "temperature", "temperature_k": 293.15})
    run_id, seconds = study.solve("TEG stack")
    picture("teg_stack_exact")
    path = save("teg_stack_exact")
    m = saved_run(path, study.id, run_id)["metrics"]
    # The heat is solved after the flow, so point goals have no value: the drop across the
    # layer is read from the copper on either side, which copper keeps nearly even, as the
    # coolest of the block above less the warmest of the block below.
    with zipfile.ZipFile(path) as z:
        document = json.loads(z.read("document.json"))
    project = document.get("multiphysics") or document.get("document", {}).get("multiphysics")
    record = next(r for r in next(s for s in project["studies"] if s["id"] == study.id)["runs"] if r["id"] == run_id)
    report = json.loads(next(a["text"] for a in record["artifacts"] if a["name"] == "thermal_report.json"))
    parts = {c["name"]: c for c in report["components"]}
    drop = parts["CopperBlock2"]["min_temperature_k"] - parts["CopperBlock"]["max_temperature_k"]
    exact = 5.0 * 0.002 / (1.5 * 0.025**2)
    return {
        "cells": study.dimensions(),
        "drop_k": drop,
        "drop_exact_k": exact,
        "parts": [{"name": n, "min_c": c["min_temperature_k"] - 273.15, "mean_c": c["mean_temperature_k"] - 273.15, "max_c": c["max_temperature_k"] - 273.15} for n, c in parts.items()],
        "error_percent": (drop / exact - 1) * 100 if exact else None,
        "total_power_w": m.get("thermal_total_power_w"),
        "relative_energy_imbalance": m.get("thermal_relative_energy_imbalance"),
        "run_seconds": seconds,
        "notes": study.notes,
    }


@case
def porous_brick():
    # A catalyst-like brick filling a 40 × 30 mm duct for 15 mm of its 150 mm, with a Darcy
    # coefficient of 1e7 1/m², air at 0.1 m/s. The flow through it is even, so the pressure falls
    # across it by Darcy's μ·D·U·L = 1.8e-5 × 1e7 × 0.1 × 0.015 = 0.27 Pa; the duct's own
    # friction over the same length adds under 0.001 Pa. Stands in for the catalyst heat-up case
    # (022), whose housing comes as an open surface and as an STL past the conversion limit.
    new_design()
    wall, w, h, length = 2.5, 40.0, 30.0, 150.0  # on the 2.5 mm cells exactly
    run("solid.box", {"length": length, "width": w + 2 * wall, "height": h + 2 * wall, "corner": [0, 0, 0], "operation": "new", "body_name": "Wall"})
    run("solid.box", {"length": length + 10, "width": w, "height": h, "corner": [-5, wall, wall], "operation": "cut"})
    run("solid.box", {"length": 15, "width": w, "height": h, "corner": [67.5, wall, wall], "operation": "new", "body_name": "Brick"})
    flow = {"kinematic_viscosity_m2_per_s": AIR_NU, "density_kg_per_m3": AIR_RHO, "end_time_s": 3.0, "converge_tolerance": 2e-4}
    study = Study("Porous brick", 2.5, flow)
    d = 1e7
    study.set(porous=[{"body": "Brick", "viscous_resistance_per_m2": d, "inertial_resistance_per_m": 0.0}])
    mid = [wall + w / 2, wall + h / 2]
    study.condition("Inlet", "Wall", [0, wall / 2, mid[1]], {"kind": "inlet", "speed_m_per_s": 0.1})
    study.condition("Outlet", "Wall", [length, wall / 2, mid[1]], {"kind": "outlet", "pressure_m2_per_s2": 0.0})
    study.goal("before", "point_pressure", point_mm=[62.5, mid[0], mid[1]])
    study.goal("after", "point_pressure", point_mm=[87.5, mid[0], mid[1]])
    study.goal("q_out", "face_volume_flow", face={"body": "Wall", "point": [length, wall / 2, mid[1]]})
    run_id, seconds = study.solve("porous brick")
    picture("porous_brick")
    m = saved_run(save("porous_brick"), study.id, run_id)["metrics"]
    q = abs(m.get("goal_q_out") or 0.0)
    u = q / (w * h * 1e-6)
    # The brick's cells: the drop is over the cells the brick holds, as many as its box covers.
    expected = AIR_NU * AIR_RHO * d * u * 0.015
    measured = (m.get("goal_before") or 0.0) - (m.get("goal_after") or 0.0)
    return {
        "cells": study.dimensions(),
        "speed_m_per_s": u,
        "drop_pa": measured,
        "drop_darcy_pa": expected,
        "error_percent": (measured / expected - 1) * 100 if expected else None,
        "goals_converged": m.get("goals_converged"),
        "simulated_s": m.get("time_s"),
        "run_seconds": seconds,
        "notes": study.notes,
    }


@case
def catalyst_flow_exact():
    # The catalyst heat-up case's exhaust duct (its walls, inlet and outlet joined into one
    # closed solid, cut out of a block) and its two 300 × 200 × 60 mm catalyst bricks, measured,
    # 200 mm apart. The bricks are porous bodies with a Darcy coefficient of 1e7 1/m² (a
    # permeability of 1e-7 m², a fine monolith); air at 2 m/s through the inlet. Heat is out of
    # reach (no heat exchange inside a porous zone), so the flow alone: the pressure drop across
    # the first brick against Darcy's μ·D·U·L, with U the flow through the brick's face.
    new_design()
    stl = os.path.join(GEOMETRY, "catalystHeatUp/cases-12/catalystHeatUp/cad")
    joined = os.path.join(OUTDIR, "catalyst_duct.stl")
    open_edges = closed_stl([os.path.join(stl, "stl", n) for n in ("walls.stl", "inlet.stl", "outlet.stl")], joined)
    print("  duct joined, edges left open:", open_edges)
    # The case's files are in metres.
    duct = run("file.insert_mesh", {"path": joined, "units": "m", "name": "Duct"})["created_bodies"][0]
    run("mesh.convert", {"body": duct, "name": "Duct solid"})
    run("solid.box", {"length": 401, "width": 210, "height": 1072, "corner": [-155, -105, -372], "operation": "new", "body_name": "Wall"})
    run("solid.combine", {"target": "Wall", "tools": [duct], "operation": "cut"})
    run("solid.box", {"length": 300, "width": 200, "height": 60, "corner": [-150, -100, 0], "operation": "new", "body_name": "SCR1"})
    run("solid.box", {"length": 300, "width": 200, "height": 60, "corner": [-150, -100, 200], "operation": "new", "body_name": "SCR2"})
    flow = {"kinematic_viscosity_m2_per_s": AIR_NU, "density_kg_per_m3": AIR_RHO, "end_time_s": 3.0, "converge_tolerance": 1e-3}
    study = Study("Catalyst flow (case geometry)", 14.0, flow)
    d = 1e7
    study.set(porous=[{"body": b, "viscous_resistance_per_m2": d, "inertial_resistance_per_m": 0.0} for b in ("SCR1", "SCR2")])
    # The inlet: the pipe's end, at x 207…241, z −372…−278 (a slanted face: a patch on the wall).
    study.condition("Inlet", "Wall", [224, 0, -372], {"kind": "inlet", "speed_m_per_s": 2.0})
    study.condition("Outlet", "Wall", [80, 102, 700], {"kind": "outlet", "pressure_m2_per_s2": 0.0})
    study.goal("below", "point_pressure", point_mm=[0, 0, -10])
    study.goal("above", "point_pressure", point_mm=[0, 0, 70])
    study.goal("q_out", "face_volume_flow", face={"body": "Wall", "point": [80, 102, 700]})
    run_id, seconds = study.solve("catalyst")
    picture("catalyst_flow_exact")
    m = saved_run(save("catalyst_flow_exact"), study.id, run_id)["metrics"]
    q = abs(m.get("goal_q_out") or 0.0)
    u = q / (0.3 * 0.2)
    expected = AIR_NU * AIR_RHO * d * u * 0.06
    measured = (m.get("goal_below") or 0.0) - (m.get("goal_above") or 0.0)
    return {
        "cells": study.dimensions(),
        "flow_m3_per_s": q,
        "brick_speed_m_per_s": u,
        "drop_pa": measured,
        "drop_darcy_pa": expected,
        "error_percent": (measured / expected - 1) * 100 if expected else None,
        "goals_converged": m.get("goals_converged"),
        "simulated_s": m.get("time_s"),
        "run_seconds": seconds,
        "notes": study.notes,
    }


print(json.dumps(RESULTS, indent=2))
