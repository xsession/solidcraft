"""Measure an STL: its connected pieces, each piece's bounding box, triangle count and volume.

Used to rebuild a training case's parts in SolidCraft from their measured sizes and places
(a battery cell as a cylinder, a channel as a box). Reads ASCII and binary STL.

    python3 stl_parts.py FILE.stl [--scale 1000]
"""
import struct
import sys

import numpy as np


def load(path):
    with open(path, "rb") as f:
        data = f.read()
    if data[:5].lower() == b"solid" and b"facet" in data[:400]:
        values = [float(t) for line in data.decode("ascii", "replace").splitlines() if line.strip().startswith("vertex") for t in line.split()[1:4]]
        return np.array(values, dtype=float).reshape(-1, 3, 3)
    count = struct.unpack("<I", data[80:84])[0]
    raw = np.frombuffer(data[84:84 + count * 50], dtype=np.uint8).reshape(count, 50)
    return raw[:, 12:48].copy().view("<f4").reshape(count, 3, 3).astype(float)


def pieces(tri):
    """Connected pieces by shared vertices (rounded to 1e-6 of the model's size)."""
    span = np.ptp(tri.reshape(-1, 3), axis=0).max() or 1.0
    keys = np.round(tri.reshape(-1, 3) / (span * 1e-6)).astype(np.int64)
    _, ids = np.unique(keys, axis=0, return_inverse=True)
    ids = ids.reshape(-1, 3)
    parent = np.arange(ids.max() + 1)

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for a, b, c in ids:
        ra, rb, rc = find(a), find(b), find(c)
        parent[rb] = ra
        parent[find(rc)] = ra
    roots = np.array([find(v) for v in ids[:, 0]])
    return [tri[roots == r] for r in np.unique(roots)]


def volume(tri):
    return float(np.einsum("ij,ij->i", tri[:, 0], np.cross(tri[:, 1], tri[:, 2])).sum() / 6.0)


if __name__ == "__main__":
    path = sys.argv[1]
    scale = float(sys.argv[sys.argv.index("--scale") + 1]) if "--scale" in sys.argv else 1.0
    tri = load(path) * scale
    lo, hi = tri.reshape(-1, 3).min(0), tri.reshape(-1, 3).max(0)
    print(f"{path}: {len(tri)} triangles, box {np.round(lo, 3)} .. {np.round(hi, 3)}, volume {volume(tri):.6g}")
    parts = pieces(tri)
    print(f"{len(parts)} pieces")
    for i, p in enumerate(sorted(parts, key=lambda p: tuple(p.reshape(-1, 3).min(0)))):
        a, b = p.reshape(-1, 3).min(0), p.reshape(-1, 3).max(0)
        print(f"  {i}: {len(p)} tri, min {np.round(a, 3)}, max {np.round(b, 3)}, size {np.round(b - a, 3)}, volume {volume(p):.6g}")
