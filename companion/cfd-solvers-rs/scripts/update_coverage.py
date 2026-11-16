#!/usr/bin/env python3
"""Rebuild the source-unit ledger; source presence never implies qualification."""
import argparse
import collections
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def coverage_markdown(report):
    groups = collections.defaultdict(collections.Counter)
    for unit in report['units']:
        groups[unit['group']][unit['status']] += 1
    out = ['# Source-unit coverage', '', f"Source revision: `{report['upstream']['revision']}`.", '',
           '**Full clone is unfinished.** Counts refer to public header units, not solvers or effort. Native reference means the scope stated below exists, with no whole-API/backend qualification.', '',
           '| Group | Native reference | Partial | Not ported |', '| --- | ---: | ---: | ---: |']
    for group, counts in sorted(groups.items()):
        out.append(f"| {group} | {counts['native_reference']} | {counts['native_partial']} | {counts['not_ported']} |")
    out += ['', '## Every public unit', '', '| Upstream public header | Status | Rust scope / remaining boundary |', '| --- | --- | --- |']
    for unit in report['units']:
        scope = unit.get('scope', 'No native implementation.').replace('|', '/')
        out.append(f"| [{unit['path']}]({unit['source_url']}) | {unit['status']} | {scope} |")
    return '\n'.join(out) + '\n'

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--check', action='store_true')
    args = parser.parse_args()
    inventory = json.loads((ROOT / 'parity/upstream-units.json').read_text())
    port_map = json.loads((ROOT / 'parity/port-map.json').read_text())
    known = {x['path'] for x in inventory['units']}
    unknown = set(port_map) - known
    if unknown:
        raise SystemExit(f'Unknown upstream units: {sorted(unknown)}')
    units = []
    for original in inventory['units']:
        entry = dict(original)
        entry.update(status='not_ported', native_files=[], limitations=['No native implementation in this release.'])
        entry.update(port_map.get(entry['path'], {}))
        entry['fully_qualified_parity'] = False
        for filename in entry['native_files']:
            if not (ROOT / filename).is_file():
                raise SystemExit(f'Missing native file: {filename}')
        units.append(entry)
    counts = dict(collections.Counter(u['status'] for u in units))
    report = {
        'schema': 'cfd-solvers-rs/coverage/1',
        'upstream': inventory['upstream'],
        'complete_clone': False,
        'meaning': 'Public header units are inventory entries, not independent solvers. Native reference implementation is not full API/backend/physical qualification.',
        'unit_count': len(units), 'status_counts': counts,
        'backends': {'serial_cpu': 'implemented_for_listed_methods', 'rayon_cpu': 'D2Q9 cell updates only', 'sycl': 'not_ported', 'mpi': 'not_ported', 'cuda': 'not_ported'},
        'units': units,
    }
    output = json.dumps(report, indent=2) + '\n'
    target = ROOT / 'parity/status.json'
    markdown = coverage_markdown(report)
    docs = ROOT / 'docs/COVERAGE.md'
    if args.check:
        if target.read_text() != output:
            raise SystemExit('Coverage ledger is stale; run scripts/update_coverage.py')
        if docs.read_text() != markdown:
            raise SystemExit('Coverage documentation is stale; run scripts/update_coverage.py')
    else:
        target.write_text(output)
        docs.write_text(markdown)
    print(json.dumps({'unit_count': len(units), 'status_counts': counts, 'complete_clone': False}))

if __name__ == '__main__':
    main()
