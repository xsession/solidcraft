"""Write example_cases/templates/README.md: every case template SolidCraft ships, the solver it
runs and what it shows, from crates/multiphysics/src/templates.rs and the template files, and,
when a results file from the all-templates run is given, whether it ran and what it reported.

    python3 make_index.py [results.json]

results.json: {template id: {"completed": bool, "seconds": float, "metrics": {...}}}, as
summarise_runs.py writes it from the all-templates runner's output folder.
"""
import json
import os
import re
import sys

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.normpath(os.path.join(here, "..", ".."))
source = open(os.path.join(root, "crates/multiphysics/src/templates.rs")).read()
entries = re.findall(r'\(\s*"([a-z0-9-]+)",\s*"([^"]+)",\s*include_str!\("\.\./data/templates/([a-z0-9-]+\.json)"\)', source)
results = json.load(open(sys.argv[1])) if len(sys.argv) > 1 else {}

rows = []
for template_id, name, file in entries:
    case = json.load(open(os.path.join(root, "crates/multiphysics/data/templates", file)))
    family, _, shows = name.partition(" · ")
    result = results.get(template_id)
    if result is None:
        ran = "–"
    elif result.get("completed"):
        ran = f"yes, {result.get('seconds', 0):.1f} s"
    else:
        ran = "**no**: " + result.get("error", "did not complete")
    rows.append((family, template_id, shows or name, case.get("solver", "?"), ran))

rows.sort()
lines = [
    "# Case templates",
    "",
    "Every case template SolidCraft offers for a new study (**New project** in the Multiphysics",
    "workspace, or `multiphysics.study.create {\"template\": …}`), each a small, complete case of",
    "one solver of the companion. `real_reference_runner_all_templates` in",
    "`apps/multiphysics-host` runs all of them through a real companion; this page is written by",
    "`make_index.py` from the templates themselves and that run's results.",
    "",
    f"{len(rows)} templates.",
    "",
    "| Family | Template | Shows | Solver | Ran |",
    "|---|---|---|---|---|",
]
lines += [f"| {f} | `{t}` | {s} | `{v}` | {r} |" for f, t, s, v, r in rows]
open(os.path.join(here, "README.md"), "w").write("\n".join(lines) + "\n")
print(f"{len(rows)} templates")
