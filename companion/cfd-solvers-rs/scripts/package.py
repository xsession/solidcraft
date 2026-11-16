#!/usr/bin/env python3
"""Package the source, vendored crates and any staged executable/sample results."""
import argparse
import hashlib
from pathlib import Path
import zipfile

ROOT = Path(__file__).resolve().parents[1]
SKIP = {'target', 'dist', '.git', '__pycache__'}

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    output = a.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    files = [f for f in sorted(ROOT.rglob('*')) if f.is_file()
             and not any(part in SKIP for part in f.relative_to(ROOT).parts)
             and f.resolve() != output
             and f.name != 'package-files.sha256']
    checksum = ROOT/'parity/package-files.sha256'
    checksum.write_text(''.join(hashlib.sha256(f.read_bytes()).hexdigest()+'  '+f.relative_to(ROOT).as_posix()+'\n' for f in files))
    files.append(checksum)
    with zipfile.ZipFile(output, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for file in files:
            archive.write(file, Path('cfd-solvers-rs')/file.relative_to(ROOT))
    with zipfile.ZipFile(output) as archive:
        bad = archive.testzip()
        if bad:
            raise SystemExit(f'Archive integrity error: {bad}')
    print(f'{output.name}: {len(files)} files, {output.stat().st_size} bytes, SHA256 {hashlib.sha256(output.read_bytes()).hexdigest()}')

if __name__ == '__main__':
    main()
