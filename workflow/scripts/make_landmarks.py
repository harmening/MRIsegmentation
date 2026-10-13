#!/usr/bin/env python3
"""Build the 10-10 landmark system on a head model's scalp surface.

    make_landmarks.py --scalp headmodel/mask_scalp.obj --reference T1.nii
        --fiducials config/landmarks/sub01.mrk.json
        --out headmodel/landmarks.mrk.json [--snap Nz Iz LPA RPA] [--max-snap-mm 20]

`--fiducials` is a 3D Slicer markup file holding `Nz`, `Iz`, `LPA` and `RPA`
(`Cz` optional - the builder re-estimates the cranial vertex regardless) in the
RAS of `--reference`, in mm. Where those coordinates come from is not this
script's business; see config/landmarks/README.md.

cedalion's LandmarksBuilder1010 walks the scalp surface between the fiducials
(:cite:`Oostenveld2001`), so the surface has to be the one the head model ships
- the same surface optodes are later placed on - and in mm with a superior z,
which is why the mesh is transformed from voxel indices to RAS here.

The fiducials are snapped onto the scalp first (`--snap`, a subset of the
labels; empty for none). A fiducial further than `--max-snap-mm` from the
surface is an error: that is what a file in the wrong space looks like.

Needs cedalion, nibabel, numpy, trimesh.
"""
import argparse
import json
import os
import sys

os.environ.setdefault("TQDM_DISABLE", "1")

import cedalion  # noqa: E402
import cedalion.dataclasses as cdc  # noqa: E402
import nibabel as nib  # noqa: E402
import numpy as np  # noqa: E402
import trimesh  # noqa: E402
import xarray as xr  # noqa: E402
from cedalion.geometry.landmarks import LandmarksBuilder1010  # noqa: E402
from cedalion.io.probe_geometry import read_mrk_json, save_mrk_json  # noqa: E402

REQUIRED = ["Nz", "Iz", "LPA", "RPA"]
# read_mrk_json ignores the file's coordinateSystem, so check it here: 3D Slicer
# writes LPS by default, which would silently mirror x and y
RAS_LIKE = {"RAS", "MNI", "MNI152", "MNI305", "SCANNER_RAS", "WORLD"}


def scalp_in_ras(obj_file, reference):
    """The scalp OBJ (1-based voxel indices) as a surface in RAS mm."""
    mesh = trimesh.load(obj_file, process=False)
    # the affine maps 0-based voxel centres to world coordinates, so undo the
    # head model's 1-based indexing first (as cedalion does when it loads them)
    ras = nib.affines.apply_affine(
        np.asarray(reference.affine, float), np.asarray(mesh.vertices) - 1.0
    )
    mesh = trimesh.Trimesh(vertices=ras, faces=mesh.faces, process=False)
    return cdc.TrimeshSurface(mesh, crs="ras", units=cedalion.units.Unit("mm"))


def read_fiducials(path):
    """The fiducials of a markup file, in RAS mm, with their stated CRS name."""
    with open(path) as fh:
        markups = json.load(fh)["markups"]
    systems = {str(m.get("coordinateSystem", "")).upper() for m in markups}
    if len(systems) != 1:
        sys.exit(f"{path}: expected one coordinateSystem, found {sorted(systems)}")
    system = systems.pop()

    points = read_mrk_json(path, crs="ras")
    if system == "LPS":
        # LPS -> RAS: the first two axes point the other way
        print(f"{path}: coordinateSystem LPS, converting to RAS")
        points = points * xr.DataArray([-1.0, -1.0, 1.0], dims=["ras"])
    elif system not in RAS_LIKE:
        sys.exit(
            f"{path}: coordinateSystem {system!r} is neither LPS nor one of "
            f"{sorted(RAS_LIKE)}; the fiducials must be in the RAS of the reference image"
        )
    return points, system


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--scalp", required=True, help="scalp surface (OBJ, 1-based voxel indices)")
    p.add_argument("--reference", required=True, help="NIfTI defining the voxel grid of the mesh")
    p.add_argument("--fiducials", required=True, help="markup file with Nz, Iz, LPA, RPA (+ Cz)")
    p.add_argument("--out", required=True)
    p.add_argument("--snap", nargs="*", default=REQUIRED,
                   help="fiducials to snap onto the scalp first (default all four, none if empty)")
    p.add_argument("--max-snap-mm", type=float, default=20.0,
                   help="a fiducial further than this from the scalp is an error")
    args = p.parse_args()

    scalp = scalp_in_ras(args.scalp, nib.load(args.reference))
    print(f"{args.scalp}: {scalp.nvertices} vertices in RAS mm")

    fiducials, system = read_fiducials(args.fiducials)
    labels = [str(x) for x in fiducials.label.values]
    missing = [f for f in REQUIRED if f not in labels]
    if missing:
        sys.exit(f"{args.fiducials} is missing the fiducial(s) {missing}")
    print(f"{args.fiducials}: {labels} ({system})")

    # distance to the scalp, both as the snapping amount and as a check that the
    # fiducials live in the same space as the surface
    snapped = scalp.snap(fiducials)
    offsets = np.linalg.norm(
        np.asarray(snapped.pint.dequantify().values) - np.asarray(fiducials.pint.dequantify().values),
        axis=1,
    )
    for label, offset in zip(labels, offsets):
        print(f"  {label:<4} {offset:6.2f} mm from the scalp" + ("  -> snapped" if label in args.snap else ""))
    if offsets.max() > args.max_snap_mm:
        worst = labels[int(np.argmax(offsets))]
        sys.exit(
            f"fiducial {worst} is {offsets.max():.1f} mm from the scalp "
            f"(--max-snap-mm {args.max_snap_mm}); are the fiducials in the RAS of "
            f"{os.path.basename(args.reference)}?"
        )

    if args.snap:
        keep = [label not in args.snap for label in labels]
        fiducials = xr.concat(
            (fiducials[keep], snapped[[not k for k in keep]]), dim="label"
        )

    landmarks = LandmarksBuilder1010(scalp, fiducials).build()
    save_mrk_json(args.out, landmarks, system.lower())
    built = [str(x) for x in landmarks.label.values]
    print(f"wrote {args.out}: {len(built)} landmarks, coordinateSystem {system.lower()}")


if __name__ == "__main__":
    main()
