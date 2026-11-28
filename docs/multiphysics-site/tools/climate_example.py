"""Set up the climate tutorial's study: an open-hardware weather station's board in its printed
case, a year of temperature, humidity and dust.

The design is the Universal Weather Node v2.0 (CERN-OHL-P-2.0): its KiCad board and the STEP
assembly of its 3D-printed case, from https://github.com/franzenhof/universal-weather-node
(pcb/v2.0/ and 3dprint/v20-v4/). Fetch those two folders, then:

    SOLIDCRAFT_CLI=/path/to/solidcraft-cli python3 climate_example.py /path/to/universal-weather-node work/

It writes work/weather-node.sc with the study "Weather node climate", ready for
climate_year.py prepare. Every step is a SolidCraft command; the comments say why.
"""
import json
import os
import subprocess
import sys

CLI = os.environ.get("SOLIDCRAFT_CLI", "solidcraft-cli")


def run(design_in, commands, design_out):
    """Run commands on a design (or on nothing) and save it; their results."""
    os.makedirs(os.path.dirname(os.path.abspath(design_out)), exist_ok=True)
    script = os.path.abspath(design_out + ".script.json")
    out = script + ".out"
    with open(script, "w") as f:
        json.dump({"commands": [{"command": c, "params": p} for c, p in commands]}, f)
    args = [CLI, "run", script, "--results", "--json-out", out, "--quiet", "--save", design_out]
    if design_in:
        args += ["--in", design_in]
    done = subprocess.run(args, capture_output=True, text=True)
    if done.returncode != 0:
        raise SystemExit(f"{CLI} failed: {done.stderr.strip()[-3000:]}")
    with open(out) as f:
        results = json.load(f)
    os.remove(script)
    os.remove(out)
    return results


def main(source, work):
    design = os.path.abspath(os.path.join(work, "weather-node.sc"))
    board = os.path.join(source, "pcb", "v2.0", "universal-weather-node.kicad_pcb")
    case = os.path.join(source, "3dprint", "v20-v4", "universal-weather-node-v2.0-v4-assembly.step")

    # 1. The board from KiCad (copper and silkscreen left out: a flow grid cannot hold them),
    #    then the case assembly, which also holds the project's simplified blocks for the radio
    #    module, the RJ12 jacks and the screw terminals.
    opened = run(None, [("doc.open", {"path": board, "copper": "none", "silkscreen": False}),
                        ("file.insert_step", {"path": case})], design)
    board_occurrence, case_occurrence = opened[0]["occurrence"], opened[1]["occurrence"]
    # 2. Put the board where the case expects it: the assembly's board block is at x −50..50,
    #    y 50..150, z 8.4 (its top face at 10). The KiCad board opens with its corner at the
    #    origin. Moving one placement moves the rigid group, so the case is put back after.
    names = {"COMPOUND": "Case", "COMPOUND (2)": "Lid", "COMPOUND (3)": "Button", "COMPOUND (4)": "PCB envelope",
             "COMPOUND (5)": "U1", "COMPOUND (6)": "J6", "COMPOUND (7)": "J4", "COMPOUND (8)": "J3", "COMPOUND (9)": "J12",
             "COMPOUND (10)": "J13", "COMPOUND (11)": "J11", "COMPOUND (12)": "J10", "COMPOUND (13)": "J9",
             "COMPOUND (14)": "J7", "COMPOUND (15)": "J5", "COMPOUND (16)": "J17", "COMPOUND (17)": "J14", "COMPOUND (18)": "J18",
             "universal-weather-node": "PCB"}
    commands = [("occurrence.ground", {"occurrence": board_occurrence, "grounded": False}),
                ("occurrence.move", {"occurrence": board_occurrence, "space": "world", "translate": [-50, 50, 8.4]}),
                ("occurrence.ground", {"occurrence": case_occurrence, "grounded": False}),
                ("occurrence.move", {"occurrence": case_occurrence, "space": "world", "translate": [50, -50, -8.4]})]
    commands += [("body.rename", {"body": a, "name": b}) for a, b in names.items()]
    # 3. The flow grid needs the board, not its 600 drilled holes and pads: a plain block where
    #    the board is stands in for it in the study (the KiCad board stays in the design).
    commands.append(("solid.box", {"length": 100, "width": 100, "height": 1.6, "corner": [-50, 50, 8.4],
                                   "operation": "new", "body_name": "Board"}))
    connectors = ["J3", "J4", "J5", "J6", "J7", "J9", "J10", "J11", "J12", "J13", "J14", "J17", "J18"]
    bodies = ["Case", "Lid", "Button", "Board", "U1"] + connectors
    # 4. A study of the air around the node (external, 8 mm of air on every side), air at rest
    #    to start, run for 30 s of the air's time, stopping once it settles.
    commands.append(("multiphysics.model.create", {"name": "Weather node climate", "analysis": "external", "bodies": bodies,
                                                   "margin_mm": 8, "flow": {"kinematic_viscosity_m2_per_s": 1.5e-5,
                                                                            "density_kg_per_m3": 1.2, "end_time_s": 30,
                                                                            "converge_tolerance": 1e-3}}))
    run(design, commands, design)
    sid = next(s["id"] for s in run(design, [("multiphysics.study.list", {})], design)[0]["project"]["studies"]
               if s["name"] == "Weather node climate")

    def edit(command, params):
        revision = next(s["revision"] for s in run(design, [("multiphysics.study.list", {})], design)[0]["project"]["studies"]
                        if s["id"] == sid)
        return run(design, [(command, dict(params, id=sid, expected_revision=revision))], design)[0]

    # 5. Heat with natural convection; the case is PETG, the board FR-4 with copper planes, the
    #    radio module (ESP32-S3 with Wi-Fi on and LoRa listening) 0.42 W, the regulator and the
    #    rest 0.08 W on the board.
    edit("multiphysics.model.heat", {"enabled": True, "buoyancy": True})
    edit("multiphysics.model.material", {"bodies": ["Case", "Lid", "Button"], "conductivity_w_per_m_k": 0.2})
    edit("multiphysics.model.material", {"bodies": ["Board"], "conductivity_w_per_m_k": 3.0, "power_w": 0.08,
                                         "max_temperature_k": 378.15})
    edit("multiphysics.model.material", {"bodies": ["U1"], "conductivity_w_per_m_k": 5.0, "power_w": 0.42,
                                         "max_temperature_k": 358.15})
    edit("multiphysics.model.material", {"bodies": connectors, "conductivity_w_per_m_k": 0.25})
    edit("multiphysics.model.goal", {"name": "U1_max", "kind": "body_max_temperature", "body": "U1"})
    edit("multiphysics.model.goal", {"name": "Board_mean", "kind": "body_mean_temperature", "body": "Board"})
    # 6. The case's walls (1.4 to 2.6 mm) are thinner than the cells: close every grid face they
    #    cross, so the case keeps only its real openings, the connector and cable cut-outs.
    edit("multiphysics.model.set", {"close_thin_walls": True})
    # 7. A shelf in a shed: a breeze along x, open above, a floor below.
    edit("multiphysics.model.surroundings", {"breeze_m_per_s": 0.1, "breeze_from": "-x", "floor": True})
    # 8. Farm dust, counted where it comes into the case: a box whose sides lie inside the walls.
    edit("multiphysics.model.dust", {"concentration_mg_per_m3": 0.5, "particles_per_size": 1000, "max_time_s": 3600,
                                     "entry": {"kind": "box", "low_mm": [-50.5, 21, 1], "high_mm": [50.5, 152, 38]}})
    # 9. A temperate year: four seasons, each a night, morning, afternoon and evening.
    edit("multiphysics.model.climate", {"seasons": [
        {"name": "winter", "low_c": -4, "high_c": 3, "dew_point_c": -5, "dust_mg_per_m3": 0.1},
        {"name": "spring", "low_c": 4, "high_c": 16, "dew_point_c": 3, "dust_mg_per_m3": 0.4},
        {"name": "summer", "low_c": 14, "high_c": 30, "dew_point_c": 13, "dust_mg_per_m3": 0.8},
        {"name": "autumn", "low_c": 6, "high_c": 16, "dew_point_c": 7, "dust_mg_per_m3": 0.5}]})
    print(f"{design}: study {sid} is ready; next: climate_year.py prepare {design} \"Weather node climate\" {work}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise SystemExit(__doc__)
    main(sys.argv[1], sys.argv[2])
