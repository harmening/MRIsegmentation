#!/usr/bin/env python3
"""Assemble a cedalion head model directory from the FreeSurfer path's output.

    make_headmodel.py --out-dir results/sub01/headmodel --reference T1.nii
        --high-res cortex_pial.obj --reduced cortex_pial_25000.obj
        --vertexinfo cortex_pial_25000_vertexinfo.npz --scalp scalp_smoothed_10000.obj
        --masks mask_csf.nii mask_gray.nii mask_skin.nii mask_bone.nii mask_white.nii
        --atlas-dir Schaefer2018_Parcellations --atlas Schaefer17_600
        --mni-from talairach --talairach transforms/talairach.xfm

The directory follows the layout cedalion's `TwoSurfaceHeadModel.from_surfaces`
and `cedalion.data.get_*_headmodel_files()` expect (as in the published
hm_colin27 / hm_icbm152 datasets):

    mask_{csf,gray,white,bone,skin}.nii  tissue masks of the MATLAB path
    cortex_pial_high.obj                 unreduced FreeSurfer cortex surface
    mask_brain.obj                       the decimated cortex surface
    cortex_pial_inflated.obj             its vertices on the inflated cortex
    mask_scalp.obj                       the decimated scalp surface
    brain_vertex_coordinates.csv         per vertex: parcel, MNI152 position and
                                         the vertex it came from on cortex_pial_high
    voxel_to_vertex_brain.mtx.gz         cortical ribbon voxel -> cortex vertex
    parcel_colors.json                   parcel name -> 8-bit RGB

landmarks.mrk.json, the remaining file of a complete model, is written into the
same directory by make_landmarks.py.

**Coordinates.** The meshes keep this pipeline's 1-based voxel indices, which is
what cedalion expects of a head model's surface files: it shifts them into its
own 0-based ijk itself (`TwoSurfaceHeadModel.from_surfaces`, see its
`surfaces_are_zero_based`). The MNI coordinates below are computed from the
0-based position, i.e. after that same shift, so that they describe the vertex
where cedalion will put it.

Needs cedalion (for its MNI305->MNI152 transform), nibabel, numpy, pandas, scipy.
"""
import argparse
import gzip
import json
import os
import shutil
import sys

import cedalion
import nibabel as nib
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse
from cedalion.dot.utils import mni305_to_mni152
from cedalion.io.anatomy import trimesh_from_freesurfer
from scipy.spatial import KDTree

# the same MNI .xfm parser the head mask cleanup uses
from fs_clean_headmask import read_xfm
from fs_parcels import BACKGROUND_PREFIX, parcel_colors

# tissue masks under the names cedalion's HeadModelFiles maps to its tissue types
MASKS = ("mask_csf.nii", "mask_gray.nii", "mask_skin.nii", "mask_bone.nii", "mask_white.nii")


def read_obj(path):
    """Vertices and 0-based faces of an OBJ file."""
    verts, faces = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("v "):
                verts.append(line.split()[1:4])
            elif line.startswith("f "):
                faces.append([tok.split("/")[0] for tok in line.split()[1:4]])
    return np.array(verts, float), np.array(faces, int) - 1


def write_obj(path, vertices, faces, comment):
    with open(path, "w") as fh:
        fh.write(f"# {comment}\n")
        for v in vertices:
            fh.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        for f in faces:
            fh.write(f"f {f[0] + 1} {f[1] + 1} {f[2] + 1}\n")


def mni152_from_voxel(reference, mode, talairach_file):
    """4x4 transform from 0-based voxel indices of the reference to MNI152 RAS (mm).

    `mode` says what the reference's own RAS already is:
      talairach  an individual scanner RAS, registered to MNI305 by recon-all
      mni305     MNI305 (e.g. Colin27, which is distributed registered to it)
      mni152     MNI152 already (e.g. the ICBM152 templates), nothing to do
    """
    vox2ras = np.asarray(reference.affine, float)
    if mode == "mni152":
        return vox2ras
    if mode == "talairach":
        if not talairach_file or not os.path.exists(talairach_file):
            sys.exit(f"--mni-from talairach needs the subject's talairach.xfm ({talairach_file})")
        # recon-all's talairach.xfm registers the subject's scanner RAS to MNI305
        vox2mni305 = read_xfm(talairach_file) @ vox2ras
    elif mode == "mni305":
        vox2mni305 = vox2ras
    else:
        sys.exit(f"unknown --mni-from {mode!r}")
    # mm -> mm, so the magnitude is the transform; .values would strip the unit
    # with a warning
    data = mni305_to_mni152.data
    to_mni152 = np.asarray(getattr(data, "magnitude", data), float)
    if to_mni152.shape != (4, 4):
        sys.exit(f"unexpected mni305_to_mni152 shape {to_mni152.shape}")
    return to_mni152 @ vox2mni305


def voxel_to_vertex_matrix(info, nvertices, reference):
    """cedalion's voxel_to_vertex_brain: one row per voxel of the whole grid."""
    v2v = np.asarray(info["voxel_to_vertex"])
    flat_index = np.asarray(info["voxel_flat_index"])
    grid_shape = tuple(int(n) for n in info["grid_shape"])
    if grid_shape != reference.shape[:3]:
        sys.exit(f"the mapping was built on a {grid_shape} grid, the reference is {reference.shape[:3]}")

    mapped = v2v >= 0
    rows, cols = flat_index[mapped], v2v[mapped]
    if len(np.unique(rows)) != len(rows):
        sys.exit("a voxel is mapped to more than one vertex")
    return scipy.sparse.coo_array(
        (np.ones(len(rows)), (rows, cols)),
        shape=(int(np.prod(grid_shape)), nvertices),
    )


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--out-dir", required=True)
    p.add_argument("--reference", required=True, help="NIfTI defining the voxel grid of the meshes")
    p.add_argument("--high-res", required=True, help="unreduced cortex surface (OBJ)")
    p.add_argument("--reduced", required=True, help="decimated cortex surface (OBJ)")
    p.add_argument("--vertexinfo", required=True, help="NPZ written by fs_decimate_mesh.py")
    p.add_argument("--scalp", required=True, help="decimated scalp surface (OBJ)")
    p.add_argument("--inflated", required=True, nargs=2, metavar=("LH", "RH"),
                   help="FreeSurfer lh.inflated and rh.inflated of the same recon-all run")
    p.add_argument("--masks", required=True, nargs="+", help=f"the tissue masks {MASKS}")
    p.add_argument("--atlas-dir", required=True)
    p.add_argument("--atlas", default="Schaefer17_600")
    p.add_argument("--mni-from", default="talairach", choices=("talairach", "mni305", "mni152"),
                   help="what the reference's own RAS is, see mni152_from_voxel")
    p.add_argument("--talairach", help="subject's talairach.xfm (for --mni-from talairach)")
    args = p.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    out = lambda name: os.path.join(args.out_dir, name)  # noqa: E731
    reference = nib.load(args.reference)

    # --- surfaces, in the 1-based voxel indices cedalion expects ---------------
    high_verts, high_faces = read_obj(args.high_res)
    red_verts, red_faces = read_obj(args.reduced)
    scalp_verts, scalp_faces = read_obj(args.scalp)
    for name, (verts, faces), src in (
        ("cortex_pial_high.obj", (high_verts, high_faces), args.high_res),
        ("mask_brain.obj", (red_verts, red_faces), args.reduced),
        ("mask_scalp.obj", (scalp_verts, scalp_faces), args.scalp),
    ):
        write_obj(out(name), verts, faces,
                  f"{os.path.basename(src)} in 1-based voxel indices of {os.path.basename(args.reference)}")
        print(f"wrote {name}: {len(verts)} vertices, {len(faces)} faces")

    # --- per-vertex table ------------------------------------------------------
    info = np.load(args.vertexinfo, allow_pickle=False)
    parcels = info["parcel"].astype(str)
    if len(parcels) != len(red_verts):
        sys.exit(f"{args.vertexinfo} has {len(parcels)} parcels, the mesh has {len(red_verts)} vertices")

    # the decimation keeps the positions of the vertices it retains, so every
    # reduced vertex is one of the high-resolution vertices; the distance is a
    # check on that rather than a tolerance
    dist, fsaverage_vertex = KDTree(high_verts).query(red_verts, workers=-1)
    if dist.max() > 1e-6:
        print(f"WARNING: reduced vertices moved up to {dist.max():.4f} voxels off the "
              "high-resolution surface; fsaverage_vertex is the nearest vertex, not the same one")

    # the inflated surface for parcel and activation plots: the same vertices as
    # the reduced mesh, but at their position on the inflated cortex. Its own
    # "inflated" CRS, so no voxel shift applies.
    inflated = trimesh_from_freesurfer(
        args.inflated[0], args.inflated[1],
        crs="inflated", units=cedalion.units(None).units, offset_hemispheres=True,
    )
    inflated_verts = np.asarray(inflated.mesh.vertices)
    if len(inflated_verts) != len(high_verts):
        sys.exit(f"the inflated surfaces have {len(inflated_verts)} vertices, "
                 f"{os.path.basename(args.high_res)} has {len(high_verts)}; "
                 "they must come from the same recon-all run")
    write_obj(out("cortex_pial_inflated.obj"), inflated_verts[fsaverage_vertex], red_faces,
              "inflated cortex, reduced to the vertices of mask_brain.obj")
    print(f"wrote cortex_pial_inflated.obj: {len(red_verts)} vertices, {len(red_faces)} faces")

    mni = nib.affines.apply_affine(
        mni152_from_voxel(reference, args.mni_from, args.talairach), red_verts - 1.0
    )
    pd.DataFrame({
        "vertex": np.arange(len(red_verts)),
        "fsaverage_vertex": fsaverage_vertex,
        "parcel": parcels,
        "mni152_r": mni[:, 0],
        "mni152_a": mni[:, 1],
        "mni152_s": mni[:, 2],
    }).to_csv(out("brain_vertex_coordinates.csv"), index=False, float_format="%.3f")
    background = np.char.startswith(parcels, BACKGROUND_PREFIX)
    print(f"wrote brain_vertex_coordinates.csv: {len(red_verts)} vertices, "
          f"{len(set(parcels[~background]))} parcels plus {background.sum()} "
          f"medial-wall vertices, MNI152 via {args.mni_from}")

    # --- voxel to vertex mapping ----------------------------------------------
    m = voxel_to_vertex_matrix(info, len(red_verts), reference)
    with gzip.open(out("voxel_to_vertex_brain.mtx.gz"), "wb") as fh:
        scipy.io.mmwrite(fh, m)
    unsupported = len(red_verts) - len(np.unique(m.col))
    print(f"wrote voxel_to_vertex_brain.mtx.gz: {m.shape[0]} voxels x {m.shape[1]} vertices, "
          f"{m.nnz} mapped, {unsupported} vertices without voxels")

    # --- parcel colours and tissue masks --------------------------------------
    colors = parcel_colors(args.atlas_dir, args.atlas)
    missing = sorted(set(parcels[~background]) - set(colors))
    if missing:
        sys.exit(f"no colour for parcel(s) {missing[:5]} in {args.atlas}")
    with open(out("parcel_colors.json"), "w") as fh:
        json.dump(colors, fh, indent=2)
    print(f"wrote parcel_colors.json: {len(colors)} parcels")

    if sorted(os.path.basename(m) for m in args.masks) != sorted(MASKS):
        sys.exit(f"--masks must be exactly {MASKS}, got {[os.path.basename(m) for m in args.masks]}")
    for src in args.masks:
        shutil.copyfile(src, out(os.path.basename(src)))
    print(f"copied {len(args.masks)} tissue masks")


if __name__ == "__main__":
    main()
