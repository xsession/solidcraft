"""Build the documentation's example design and run it for real.

A board in an air duct: an open-ended duct, a board, three heat-producing parts, two small
passives and a connector. The script drives a running SolidCraft through its control channel
(commands only, no clicks), runs the study twice with the native companion (1.0 and 0.5 m/s at
the inlet), keeps both runs and saves the design. The documentation's pictures are made from
the saved design (see make_images.md).

    SOLIDCRAFT_CFD_SOLVE=/path/to/cfd-solve solidcraft --control 7979 &
    SOLIDCRAFT_CONTROL_TOKEN=$(cat ~/.local/share/solidcraft/control/7979.token) \\
        python3 make_study.py ../docs/examples/board-in-duct.sc
"""
import json
import os
import socket
import sys
import time

PORT = int(os.environ.get("PORT", "7979"))
TOKEN = os.environ["SOLIDCRAFT_CONTROL_TOKEN"]
OUT = os.path.abspath(sys.argv[1])

sock = socket.create_connection(("127.0.0.1", PORT), timeout=600)
reader = sock.makefile("r")
counter = 0


def call(method, params=None):
    global counter
    counter += 1
    sock.sendall((json.dumps({"id": counter, "method": method, "params": params or {}, "token": TOKEN}) + "\n").encode())
    reply = json.loads(reader.readline())
    if not reply.get("ok", False):
        raise SystemExit(f"{method} failed: {json.dumps(reply)[:600]}")
    return reply["result"]


def run(command, params=None):
    return call("engine.execute", {"command": command, "params": params or {}})


call("ui.home", {"open": False})
# The duct: 160 × 60 × 40 mm, 4 mm walls, open at both ends.
run("solid.box", {"length": 160, "width": 60, "height": 40, "corner": [0, 0, 0], "operation": "new", "body_name": "Duct"})
run("solid.box", {"length": 170, "width": 52, "height": 32, "center": [80, 30, 20], "operation": "cut"})
# The board and its parts, named by reference designators.
run("solid.box", {"length": 120, "width": 44, "height": 4, "corner": [20, 8, 6], "operation": "new", "body_name": "PCB1"})
parts = [
    ("U1", [45, 30], [14, 14, 5]),
    ("U2", [80, 20], [10, 10, 5]),
    ("Q1", [105, 38], [8, 8, 5]),
    ("J1", [130, 30], [8, 20, 10]),
    ("C1", [62, 40], [3, 1.6, 1.6]),
    ("C2", [92, 40], [3, 1.6, 1.6]),
    ("R1", [62, 18], [3, 1.6, 0.8]),
]
for name, (x, y), (l, w, h) in parts:
    run("solid.box", {"length": l, "width": w, "height": h, "center": [x, y, 10 + h / 2], "operation": "new", "body_name": name})

# Tag the parts by their designators (U: ic, Q: transistor, C: capacitor…).
grouped = run("multiphysics.parts.group", {})
print("grouped", [(g["body"], g["class"]) for g in grouped["grouped"]])
run("multiphysics.parts.selection", {"name": "Heat sources", "members": {"kind": "query", "rules": [{"rule": "tag", "tag": "ic"}]}})

# A study of the model: air inside the duct, 4 mm cells (a run stays under the 16 MiB result limit).
study = run("multiphysics.model.create", {"name": "Board in duct", "analysis": "internal", "cell_mm": 4.0})["study"]
sid = study["id"]


def rev():
    return next(s for s in run("multiphysics.study.list")["project"]["studies"] if s["id"] == sid)["revision"]


face = lambda p: {"body": "Duct", "point": p}
run("multiphysics.model.condition", {"id": sid, "expected_revision": rev(), "name": "Inlet", "face": face([0, 30, 2]),
                                     "role": {"kind": "inlet", "speed_m_per_s": 1.0}})
run("multiphysics.model.condition", {"id": sid, "expected_revision": rev(), "name": "Outlet", "face": face([160, 30, 2]),
                                     "role": {"kind": "outlet", "pressure_m2_per_s2": 0.0}})
run("multiphysics.model.heat", {"id": sid, "expected_revision": rev(), "enabled": True, "inflow_temperature_k": 298.15})
for bodies, k, power, limit in [
    (["U1"], 150.0, 2.0, 358.15),
    (["U2"], 150.0, 1.2, 358.15),
    (["Q1"], 150.0, 1.0, 348.15),
    (["J1"], 0.3, 0.0, 358.15),
    (["PCB1"], 10.0, 0.0, None),
    (["C1", "C2", "R1"], 1.0, 0.0, 378.15),
]:
    params = {"id": sid, "expected_revision": rev(), "bodies": bodies, "conductivity_w_per_m_k": k, "power_w": power}
    if limit:
        params["max_temperature_k"] = limit
    run("multiphysics.model.material", params)
for name, kind in [("hottest_solid", "global_max_solid_temperature"), ("max_speed", "global_max_speed")]:
    try:
        run("multiphysics.model.goal", {"id": sid, "expected_revision": rev(), "name": name, "kind": kind})
    except SystemExit as e:
        print("goal skipped:", e)

# The passives are smaller than a cell: the solver refuses a part that covers no cell centre.
# Leave the parts under 4 mm out of the study (they keep their materials for when they are back).
run("multiphysics.model.leave_out_small", {"id": sid, "expected_revision": rev(), "below_mm": 4.0})


def solve(label):
    job = run("multiphysics.run.start", {"id": sid})
    job_id = job.get("job_id") or job.get("job", {}).get("id") or job.get("id")
    print("started", label, job_id)
    for _ in range(600):
        time.sleep(2)
        status = run("multiphysics.run.status")
        jobs = status.get("jobs", status if isinstance(status, list) else [])
        mine = next((j for j in jobs if j.get("id") == job_id), None)
        state = (mine or {}).get("state")
        if state in ("completed", "succeeded", "done", "finished") or (mine or {}).get("evidence"):
            break
        if state in ("failed", "cancelled", "error"):
            raise SystemExit(f"run {label} {state}: {json.dumps(mine)[:800]}")
    kept = run("multiphysics.run.keep", {"job_id": job_id})
    print("kept", label, [r.get("id") for r in kept.get("runs", [])][-1:])


solve("1.0 m/s")
run("multiphysics.model.condition", {"id": sid, "expected_revision": rev(), "name": "Inlet", "face": face([0, 30, 2]),
                                     "role": {"kind": "inlet", "speed_m_per_s": 0.5}})
solve("0.5 m/s")
run("file.save_as", {"path": OUT})
print("saved", OUT)
