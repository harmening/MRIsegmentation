#!/usr/bin/env python3
"""Schaefer parcels of the fsaverage6 atlas, transferred onto a subject's surface.

    fs_parcels.py --fs-subject <SUBJECTS_DIR>/<subj> --atlas-dir Schaefer2018_Parcellations
                  --atlas Schaefer17_600 --out-lh lh.annot --out-rh rh.annot

The atlas ships as an fsaverage6 label.gii per hemisphere plus a labelmap.json
holding the parcel names - already in cedalion's spelling, e.g.
`VisCent_ExStr_1_LH` and `Background+FreeSurfer_Defined_Medial_Wall_LH` - and
their colours. Parcels are transferred by nearest neighbour on the registered
sphere, which is what mri_surf2surf would do, so no precomputed per-subject
annotation is required.

Writing them as FreeSurfer .annot files makes them readable by everything that
expects recon-all's own annotations: the parcel-aware decimation
(fs_decimate_mesh.py), the head model export (make_headmodel.py) and
FreeSurfer's own tools.

Also used as a module, for the parcel names per vertex (transfer_parcels) and
the parcel colours (parcel_colors).
"""
import argparse
import json

import nibabel as nib
import numpy as np
from scipy.spatial import KDTree

# GIFTI spells the hemispheres L/R, FreeSurfer lh/rh
HEMIS = (("L", "lh"), ("R", "rh"))

# the medial wall carries no signal; cedalion skips it when mapping voxels and
# its head models ship no colour for it
BACKGROUND_PREFIX = "Background"


def load_labelmap(atlas_dir, atlas, hemi):
    """Parcel table of one hemisphere: list index -> {"name", "color"}."""
    with open(f"{atlas_dir}/fsaverage6.{hemi}.{atlas}.labelmap.json") as fh:
        return json.load(fh)


def transfer_labels(fs_subject, atlas_dir, atlas, hemi, fs_hemi, nvertices=None):
    """Atlas label index per vertex of one subject hemisphere."""
    sphere_avg = nib.load(f"{atlas_dir}/fsaverage6.{hemi}.sphere.reg.surf.gii").darrays[0].data
    labels_avg = nib.load(f"{atlas_dir}/fsaverage6.{hemi}.{atlas}.label.gii").darrays[0].data
    sphere_sub, _ = nib.freesurfer.read_geometry(f"{fs_subject}/surf/{fs_hemi}.sphere.reg")
    if nvertices is not None and len(sphere_sub) != nvertices:
        raise SystemExit(
            f"{fs_hemi}.sphere.reg has {len(sphere_sub)} vertices, surface has {nvertices}"
        )
    _, nn = KDTree(sphere_avg).query(sphere_sub, workers=-1)
    return labels_avg[nn]


def transfer_parcels(fs_subject, atlas_dir, atlas, counts=(None, None)):
    """Parcel name per vertex of the lh-then-rh surface, from fsaverage6."""
    names = []
    for (hemi, fs_hemi), nverts in zip(HEMIS, counts):
        lut = load_labelmap(atlas_dir, atlas, hemi)
        labels = transfer_labels(fs_subject, atlas_dir, atlas, hemi, fs_hemi, nverts)
        names.append(np.asarray([lut[i]["name"] for i in labels]))
    return np.concatenate(names)


def parcel_colors(atlas_dir, atlas, skip_background=True):
    """Parcel name -> 8-bit RGB, from the atlas colour table."""
    colors = {}
    for hemi, _ in HEMIS:
        for entry in load_labelmap(atlas_dir, atlas, hemi):
            if skip_background and entry["name"].startswith(BACKGROUND_PREFIX):
                continue
            colors[entry["name"]] = [int(round(255 * c)) for c in entry["color"]]
    return colors


def write_annot_file(path, labels, lut):
    """Write one hemisphere's labels as a FreeSurfer annotation."""
    ctab = np.array(
        [[*(int(round(255 * c)) for c in e["color"]), 0] for e in lut], dtype=int
    )
    names = [e["name"] for e in lut]
    nib.freesurfer.write_annot(path, labels, ctab, names, fill_ctab=True)

    # everything downstream keys on these names - the parcel-aware voxel mapping
    # skips the medial wall by name, and they end up in the head model's
    # brain_vertex_coordinates.csv - so make sure the file reads back as written
    back_labels, _, back_names = nib.freesurfer.read_annot(path, orig_ids=False)
    if [n.decode("utf8") for n in back_names] != names:
        raise SystemExit(f"{path}: parcel names do not survive the round trip")
    if not np.array_equal(back_labels, labels):
        raise SystemExit(f"{path}: vertex labels do not survive the round trip")


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--fs-subject", required=True, help="<SUBJECTS_DIR>/<subject>")
    p.add_argument("--atlas-dir", required=True, help="directory with the fsaverage6 assets")
    p.add_argument("--atlas", default="Schaefer17_600")
    p.add_argument("--out-lh", required=True)
    p.add_argument("--out-rh", required=True)
    args = p.parse_args()

    for (hemi, fs_hemi), out in zip(HEMIS, (args.out_lh, args.out_rh)):
        lut = load_labelmap(args.atlas_dir, args.atlas, hemi)
        labels = transfer_labels(args.fs_subject, args.atlas_dir, args.atlas, hemi, fs_hemi)
        write_annot_file(out, labels, lut)
        print(
            f"wrote {out}: {len(labels)} vertices, "
            f"{len(np.unique(labels))} of {len(lut)} parcels used"
        )


if __name__ == "__main__":
    main()
