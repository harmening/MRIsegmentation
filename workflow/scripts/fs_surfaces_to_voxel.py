#!/usr/bin/env python3
"""Export FreeSurfer surfaces as OBJ files in the voxel grid of a NIfTI image.

    fs_surfaces_to_voxel.py <image.nii> <T1.mgz> <out_dir> NAME=SURFACE ...

FreeSurfer surfaces are stored in "tkr-RAS" coordinates of the conformed volume
(T1.mgz / orig.mgz). They are mapped to scanner RAS with the conformed volume's
vox2ras * inv(vox2ras_tkr), then into voxel indices of <image.nii> with the
inverse of its affine. Vertices are written as 1-BASED voxel indices, i.e. the
convention of the MATLAB pipeline's exports (FieldTrip's mri.transform and
SPM's V.mat map 1-based indices to world coordinates, and mask_scalp.obj /
mask_brain.obj are written with their inverses).

Each NAME=SURFACE argument writes <out_dir>/NAME.obj. Faces are 1-based, as the
OBJ format requires. Combined cortex files are written when both hemispheres of
a surface type are given (lh.white + rh.white -> cortex_white.obj, LH then RH).

Needs numpy and nibabel only.
"""
import os
import sys

import nibabel as nib
import numpy as np


def tkr_to_voxel1(coords_tkr, conformed, image):
    """tkr-RAS of the conformed volume -> 1-based voxel indices of `image`."""
    hdr = conformed.header
    tkr_to_world = hdr.get_vox2ras() @ np.linalg.inv(hdr.get_vox2ras_tkr())
    world_to_vox = np.linalg.inv(image.affine)
    m = world_to_vox @ tkr_to_world
    xyz1 = np.c_[coords_tkr, np.ones(len(coords_tkr))]
    return (xyz1 @ m.T)[:, :3] + 1.0


def write_obj(path, vertices, faces):
    with open(path, "w") as fh:
        for v in vertices:
            fh.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        for f in faces:
            fh.write(f"f {f[0] + 1} {f[1] + 1} {f[2] + 1}\n")


def main(argv):
    if len(argv) < 5:
        sys.exit(__doc__)
    image = nib.load(argv[1])
    conformed = nib.load(argv[2])
    out_dir = argv[3]
    os.makedirs(out_dir, exist_ok=True)

    meshes = {}
    for spec in argv[4:]:
        name, path = spec.split("=", 1)
        coords, faces = nib.freesurfer.read_geometry(path)
        vox = tkr_to_voxel1(coords, conformed, image)
        meshes[name] = (vox, faces)
        write_obj(os.path.join(out_dir, f"{name}.obj"), vox, faces)
        print(f"{name}: {len(vox)} vertices, {len(faces)} faces")

    # combined hemispheres (LH then RH), faces of RH offset by the LH vertex count
    for kind in ("white", "pial"):
        lh, rh = meshes.get(f"lh.{kind}"), meshes.get(f"rh.{kind}")
        if lh is None or rh is None:
            continue
        vox = np.vstack([lh[0], rh[0]])
        faces = np.vstack([lh[1], rh[1] + len(lh[0])])
        write_obj(os.path.join(out_dir, f"cortex_{kind}.obj"), vox, faces)
        print(f"cortex_{kind}: {len(vox)} vertices, {len(faces)} faces")


if __name__ == "__main__":
    main(sys.argv)
