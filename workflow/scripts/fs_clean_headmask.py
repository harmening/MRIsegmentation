#!/usr/bin/env python3
"""Remove artefacts from a FreeSurfer head mask (mri_seghead output).

    fs_clean_headmask.py <seghead.mgz> <out.mgz> [options]

mri_seghead thresholds the T1 and fills, so anything bright next to the head
ends up in the mask. Two independent criteria clean that up:

  --radius R      morphological opening with a ball of R voxels, keeping the
                  largest component and restoring the eroded rim, which erases
                  thin sheets (ghosting smears attached to the scalp, as in the
                  ICBM152 template). 0 disables it.
  --talairach X   population prior: with the subject's talairach.xfm, mask
                  voxels are dropped when they lie above --min-z mm in MNI
                  space, more than --margin mm outside the MNI305 population
                  head, AND are darker than --dim-below in the (intensity
                  normalised) T1. This removes structures above the vertex that
                  touch the head and therefore survive any opening - e.g. the
                  bar across the top of the Colin27 template - while a head
                  that is simply larger than the template keeps its bright
                  scalp. Restricted to the superior part because the population
                  template has no neck.

Values are kept as written by mri_seghead. Needs numpy, scipy, nibabel
(all in fspython).
"""
import argparse
import os
import re
import sys

import nibabel as nib
import numpy as np
from scipy import ndimage


def ball(r):
    g = np.indices((2 * r + 1,) * 3) - r
    return (g ** 2).sum(0) <= r * r


def largest_component(mask):
    labels, n = ndimage.label(mask, np.ones((3, 3, 3)))
    if n == 0:
        sys.exit("empty head mask")
    return labels == (np.argmax(np.bincount(labels.ravel())[1:]) + 1)


def read_xfm(path):
    """4x4 matrix of a MNI .xfm (scanner RAS -> MNI305 RAS)."""
    text = open(path).read()
    rows = [r for r in re.search(r"Linear_Transform\s*=\s*(.*?);", text, re.S).group(1).split("\n") if r.strip()]
    m = np.vstack([np.array([float(x) for x in r.split()]) for r in rows])
    return np.vstack([m, [0, 0, 0, 1]])


def outside_population_head(mask, t1, vox2ras, xfm_file, margin, min_z, dim_below, template):
    """Mask voxels outside the MNI305 head, superior, and dark in the T1."""
    tpl = nib.load(template)
    tpl_data = np.asarray(tpl.dataobj)
    tpl_ras = tpl.header.get_vox2ras()
    env = largest_component(ndimage.binary_fill_holes(tpl_data > 20))
    # distance in mm outside the population head (template is 1 mm isotropic)
    dist = ndimage.distance_transform_edt(~env)

    ijk = np.argwhere(mask)
    ras = vox2ras @ np.c_[ijk, np.ones(len(ijk))].T
    mni = read_xfm(xfm_file) @ ras
    tvox = np.linalg.inv(tpl_ras) @ mni
    tvox = np.clip(np.round(tvox[:3]).astype(int).T, 0, np.array(tpl_data.shape) - 1)

    drop = (mni[2] > min_z) & (dist[tvox[:, 0], tvox[:, 1], tvox[:, 2]] > margin)
    if dim_below is not None:
        drop &= t1[ijk[:, 0], ijk[:, 1], ijk[:, 2]] < dim_below
    out = np.zeros_like(mask)
    out[tuple(ijk[drop].T)] = True
    return out


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("infile")
    p.add_argument("outfile")
    p.add_argument("--radius", type=int, default=1, help="opening radius in voxels (0 = off)")
    p.add_argument("--talairach", help="subject's talairach.xfm (enables the population prior)")
    p.add_argument("--t1", help="intensity normalised T1.mgz on the same grid (for --dim-below)")
    p.add_argument("--margin", type=float, default=2.0, help="mm outside the population head that is still kept")
    p.add_argument("--min-z", type=float, default=40.0, help="only apply the prior above this MNI z (mm)")
    p.add_argument("--dim-below", type=float, default=45.0,
                   help="only drop voxels darker than this in the T1 (normalised, white matter = 110); "
                        "the outer scalp rim fades to ~40-70, so higher values shave the scalp")
    p.add_argument(
        "--template",
        default=os.path.join(os.environ.get("FREESURFER_HOME", "/usr/local/freesurfer"), "average", "mni305.cor.mgz"),
    )
    args = p.parse_args()

    img = nib.load(args.infile)
    data = np.asarray(img.dataobj)
    mask = data > 0
    n0 = int(mask.sum())

    if args.talairach and args.margin >= 0:
        if not os.path.exists(args.talairach):
            print(f"WARNING: {args.talairach} not found, skipping the population prior")
        elif not os.path.exists(args.template):
            print(f"WARNING: {args.template} not found, skipping the population prior")
        else:
            t1 = np.asarray(nib.load(args.t1).dataobj) if args.t1 else None
            drop = outside_population_head(
                mask, t1, img.header.get_vox2ras(), args.talairach, args.margin, args.min_z,
                args.dim_below if t1 is not None else None, args.template,
            )
            mask &= ~drop
            print(f"head mask: {int(drop.sum())} voxels more than {args.margin:g} mm outside the MNI305 head "
                  f"above z={args.min_z:g} mm" + (f" and darker than {args.dim_below:g}" if t1 is not None else ""))

    if args.radius > 0:
        se = ball(args.radius)
        clean = mask & ndimage.binary_dilation(largest_component(ndimage.binary_opening(mask, se)), se)
    else:
        clean = largest_component(mask)
    print(f"head mask: {n0} voxels, removed {n0 - int(clean.sum())} (opening radius {args.radius})")

    nib.save(nib.MGHImage(np.where(clean, data, 0).astype(data.dtype), img.affine, img.header), args.outfile)


if __name__ == "__main__":
    main()
