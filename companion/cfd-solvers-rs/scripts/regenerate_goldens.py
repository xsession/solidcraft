#!/usr/bin/env python3
"""Compile a test-only C++ oracle against a verified upstream checkout.

Normal cargo builds and tests require no C++ compiler or checkout.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[1]
SOURCES = ['src/fem/poisson1d.cpp', 'src/lbm/d2q9.cpp',
           'src/fdtd/maxwell1d.cpp', 'src/fvm/compressible1d.cpp',
           'src/fvm/diffusion2d.cpp', 'src/fem/poisson2d.cpp', 'src/fem/poisson3d.cpp',
           'src/fem/mesh2d.cpp', 'src/fem/mesh3d.cpp', 'src/fem/reference_element.cpp', 'src/fem/heat2d.cpp',
           'src/rf/network.cpp', 'src/optics/ray.cpp', 'src/optics/sequential.cpp',
           'src/optics/gaussian_beam.cpp', 'src/optics/materials.cpp', 'src/optics/polarization.cpp',
           'src/chemistry/kinetics.cpp', 'src/chemistry/thermo.cpp', 'src/chemistry/equilibrium.cpp',
           'src/chemistry/implicit_reactor.cpp', 'src/acoustics/kspace2d.cpp', 'src/spectral/fft.cpp',
           'src/fem/elasticity2d.cpp', 'src/fem/modal_bar1d.cpp', 'src/fem/piezoelectric_bar1d.cpp', 'src/fem/material_laws.cpp']
FVM_SOURCES = ['src/core/complex_sparse.cpp', 'src/core/geometric_multigrid.cpp',
               'src/fvm/poly_mesh.cpp', 'src/fvm/operators.cpp', 'src/fvm/pressure_velocity.cpp',
               'src/fvm/schemes.cpp', 'src/fvm/runtime_schemes.cpp', 'src/fvm/scalar_transport.cpp',
               'src/fvm/collocated_incompressible.cpp', 'src/fvm/thermal_transport.cpp', 'src/fvm/solid_heat.cpp',
               'src/fvm/advanced_models.cpp', 'src/fvm/turbulence.cpp', 'src/fvm/radiation.cpp',
               'src/fvm/radiation_network.cpp', 'src/fvm/derived_fields.cpp', 'src/fvm/rans_transport.cpp']

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--upstream', type=Path, required=True)
    p.add_argument('--compiler', default=os.environ.get('CXX', 'g++'))
    p.add_argument('--check', action='store_true', help='Compare regenerated numbers without changing fixtures')
    p.add_argument('--group', choices=['baseline','fvm'], default='baseline')
    a = p.parse_args()
    upstream = a.upstream.resolve()
    sources = FVM_SOURCES if a.group == 'fvm' else SOURCES
    oracle = ROOT / ('parity/oracle_fvm.cpp' if a.group == 'fvm' else 'parity/oracle.cpp')
    suffix = '-fvm' if a.group == 'fvm' else ''
    expected = json.loads((ROOT / 'parity/upstream-files.json').read_text())
    inputs = [f for f in expected['files'] if f['path'].startswith('include/') or f['path'] in sources]
    for item in inputs:
        path = upstream / item['path']
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != item['sha256']:
            raise SystemExit(f'Upstream input does not match pinned source: {item["path"]}')
    with tempfile.TemporaryDirectory(prefix='cfd-rust-oracle-') as temp:
        exe = Path(temp) / ('oracle.exe' if os.name == 'nt' else 'oracle')
        command = [a.compiler, '-std=c++20', '-O2', '-I', str(upstream / 'include'),
                   str(oracle), *[str(upstream / s) for s in sources], '-o', str(exe)]
        subprocess.run(command, check=True)
        result = subprocess.run([str(exe)], check=True, capture_output=True, text=True)
        values = json.loads(result.stdout)
    target = ROOT / f'parity/goldens{suffix}.json'
    if a.check:
        stored = json.loads(target.read_text())
        if set(stored) != set(values):
            raise SystemExit('Golden scenario names differ')
        for name, actual in values.items():
            if len(actual) != len(stored[name]):
                raise SystemExit(f'Golden shape mismatch: {name}')
            for index, (x, y) in enumerate(zip(actual, stored[name])):
                tolerance = (2e-6 if name.startswith('lbm_') else 1e-10) * max(1.0, abs(y))
                if abs(x-y) > tolerance:
                    raise SystemExit(f'Golden numerical mismatch: {name}[{index}]: {x} vs {y}')
        print(json.dumps({'checked_scenarios': len(values), 'upstream_revision': expected['upstream']['revision']}))
        return
    target.write_text(json.dumps(values, indent=2) + '\n')
    metadata = {'revision': expected['upstream']['revision'], 'compiler': subprocess.check_output([a.compiler, '--version'], text=True).splitlines()[0],
                'cases': len(values), 'numeric_values': sum(len(v) for v in values.values()),
                'fixture_sha256': hashlib.sha256(target.read_bytes()).hexdigest(),
                'oracle_sha256': hashlib.sha256(oracle.read_bytes()).hexdigest(),
                'method': 'C++20 -O2 serial CPU; original sources compiled separately; no C++ code linked into Rust.'}
    (ROOT / f'parity/goldens{suffix}-metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')
    print(json.dumps(metadata))

if __name__ == '__main__':
    main()
