"""Summarise the all-templates run's output folder (SOLIDCRAFT_TEST_OUTPUT_DIR of
real_reference_runner_all_templates: one folder per template with its manifest) into the
results file make_index.py reads.

    python3 summarise_runs.py OUTPUT_DIR results.json
"""
import json
import os
import sys

out, target = sys.argv[1], sys.argv[2]
results = {}
for name in sorted(os.listdir(out)):
    manifest = os.path.join(out, name, "manifest.json")
    failed = os.path.join(out, name, "failed.txt")
    if os.path.exists(manifest):
        m = json.load(open(manifest))
        metrics = m.get("metrics", {})
        results[name] = {"completed": bool(m.get("execution_completed")), "seconds": metrics.get("wall_time_s", 0.0), "metrics": metrics}
    elif os.path.exists(failed):
        results[name] = {"completed": False, "error": open(failed).read().strip()[:200]}
json.dump(results, open(target, "w"), indent=1)
print(len(results), "templates;", sum(1 for r in results.values() if r.get("completed")), "completed")
