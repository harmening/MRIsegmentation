#!/usr/bin/env python3
"""QC: compare MATLAB-pipeline and FreeSurfer surfaces of one subject.

    compare_surfaces.py <image.nii> <out.png> <out.json> NAME=GROUP:COLOR:path.obj ...

All OBJ files must be in the same coordinate system (1-based voxel indices of
<image.nii>, as written by the MATLAB pipeline and fs_surfaces_to_voxel.py).
GROUP is "scalp" or "cortex"; surfaces of the same group are compared with
each other (nearest-vertex distances, in voxels), and every surface is drawn
as its exact cross-section (triangle/plane intersections, not vertex slabs)
on a sagittal, coronal and axial slice through the volume centre.

Needs numpy, scipy, nibabel, matplotlib.
"""
import json
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import nibabel as nib  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.collections import LineCollection  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402


def read_obj(path):
    verts, faces = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("v "):
                verts.append(line.split()[1:4])
            elif line.startswith("f "):
                faces.append([tok.split("/")[0] for tok in line.split()[1:4]])
    return np.array(verts, float) - 1.0, np.array(faces, int) - 1  # 0-based voxel coords


def section(verts, faces, axis, value):
    """Line segments of the mesh cut by the plane coord[axis] == value (vectorised)."""
    tri = verts[faces]  # (n, 3, 3)
    s = tri[:, :, axis] - value
    keep = [k for k in range(3) if k != axis]
    segs = []
    for i, j in ((0, 1), (1, 2), (2, 0)):
        cross = (s[:, i] < 0) != (s[:, j] < 0)
        t = s[cross, i] / (s[cross, i] - s[cross, j])
        p = tri[cross, i] + t[:, None] * (tri[cross, j] - tri[cross, i])
        segs.append((np.flatnonzero(cross), p[:, keep]))
    # collect the two crossing points of every cut triangle
    pts = {}
    for idx, p in segs:
        for k, q in zip(idx, p):
            pts.setdefault(k, []).append(q)
    return np.array([v for v in pts.values() if len(v) == 2])


def main(argv):
    if len(argv) < 5:
        sys.exit(__doc__)
    image = nib.load(argv[1])
    out_png, out_json = argv[2], argv[3]
    surfaces = {}
    for spec in argv[4:]:
        name, rest = spec.split("=", 1)
        group, color, path = rest.split(":", 2)
        surfaces[name] = dict(group=group, color=color, mesh=read_obj(path))

    # nearest-vertex distances within each group (both directions)
    stats = {}
    names = list(surfaces)
    for a in names:
        for b in names:
            if a >= b or surfaces[a]["group"] != surfaces[b]["group"]:
                continue
            va, vb = surfaces[a]["mesh"][0], surfaces[b]["mesh"][0]
            for src, dst, vs, vd in ((a, b, va, vb), (b, a, vb, va)):
                d, _ = cKDTree(vd).query(vs)
                stats[f"{src} -> {dst}"] = dict(
                    median=float(np.median(d)), p95=float(np.percentile(d, 95)), max=float(d.max()), n=int(len(d))
                )
    with open(out_json, "w") as fh:
        json.dump(stats, fh, indent=2)

    data = image.get_fdata()
    centre = [n // 2 for n in data.shape]
    fig, axes = plt.subplots(1, 3, figsize=(18, 6.5))
    for ax, (axis, title) in zip(axes, enumerate(("sagittal", "coronal", "axial"))):
        sl = np.take(data, centre[axis], axis=axis).T
        ax.imshow(sl, cmap="gray", origin="lower", vmax=np.percentile(sl, 99.5))
        for name, s in surfaces.items():
            segs = section(*s["mesh"], axis, centre[axis])
            if len(segs):
                ax.add_collection(LineCollection(segs, colors=s["color"], linewidths=1.0, label=name))
        ax.set_title(f"{title} (index {centre[axis]})")
        ax.set_xlim(0, sl.shape[1])
        ax.set_ylim(0, sl.shape[0])
    axes[0].legend(loc="lower left", fontsize=8)
    lines = [f"{k}: median {v['median']:.2f}, 95% {v['p95']:.2f} voxels" for k, v in stats.items()]
    fig.suptitle("Exact mesh cross-sections in the NIfTI voxel grid", fontsize=11)
    fig.text(0.01, 0.01, "nearest-vertex distances  " + "\n".join(lines), fontsize=7, family="monospace", va="bottom")
    plt.tight_layout(rect=(0, 0.16, 1, 0.97))
    plt.savefig(out_png, dpi=90)
    print("\n".join(lines))


if __name__ == "__main__":
    main(sys.argv)
