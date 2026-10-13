#!/usr/bin/env python3
"""QC: compare the MATLAB tissue masks with FreeSurfer's segmentation.

    compare_masks.py --reference T1.nii --masks <dir> --aseg aseg.mgz --head seghead.mgz
                     --out-png masks.png --out-json masks.json

The two pipelines segment different things: the MATLAB path (SPM + Andy's
tools) labels the whole head in six tissue classes for FEM/BEM modelling, while
recon-all labels brain anatomy (aseg) and has no skull, scalp or subarachnoid
CSF class at all - its head mask exists only because the workflow runs
mkheadsurf. This reports what can be compared:

* volumes and Dice of white matter, grey matter, CSF, brain and head,
* for every FreeSurfer class, how the MATLAB masks label those voxels, which
  is where the differences become interpretable (e.g. FreeSurfer's filled
  white matter contains tissue that SPM calls grey),
* slices of both segmentations and of their brain-mask agreement.

FreeSurfer volumes are resampled (nearest neighbour) onto the grid of
--reference, the grid the MATLAB masks live on.

Needs numpy, nibabel, matplotlib.
"""
import argparse
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import nibabel as nib  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import ListedColormap  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

# FreeSurfer aseg labels, see FreeSurferColorLUT.txt
ASEG = {
    "white matter": [2, 41, 7, 46, 251, 252, 253, 254, 255],   # cerebral + cerebellar WM, corpus callosum
    "cortical GM": [3, 42],
    "subcortical GM": [10, 11, 12, 13, 17, 18, 26, 28, 49, 50, 51, 52, 53, 54, 58, 60],
    "cerebellar GM": [8, 47],
    "brainstem": [16],
    "ventricles/CSF": [4, 5, 14, 15, 24, 43, 44, 72],
}
MATLAB_TISSUES = ["white", "gray", "csf", "bone", "skin"]
# label value -> (name, colour) for the overview figure
PALETTE = [(0, 0, 0, 0), (1, .85, .6), (.6, .6, .6), (.3, .6, 1), (.9, .3, .3), (1, .95, .4)]
PALETTE_NAMES = ["scalp", "skull", "CSF", "grey", "white"]


def resample(volume_file, reference):
    """Nearest-neighbour resampling onto the reference grid."""
    vol = nib.load(volume_file)
    data = np.asarray(vol.dataobj)
    m = np.linalg.inv(vol.affine) @ reference.affine
    grid = np.indices(reference.shape[:3]).reshape(3, -1)
    idx = np.round(m[:3, :3] @ grid + m[:3, 3:]).astype(int)
    ok = np.all((idx >= 0) & (idx < np.array(data.shape[:3])[:, None]), axis=0)
    out = np.zeros(grid.shape[1], data.dtype)
    out[ok] = data[idx[0, ok], idx[1, ok], idx[2, ok]]
    return out.reshape(reference.shape[:3])


def dice(a, b):
    total = int(a.sum()) + int(b.sum())
    return 2.0 * int((a & b).sum()) / total if total else float("nan")


def overview_label_maps(mat, fs_class):
    """Winner-takes-all label volumes for the figure (see PALETTE)."""
    m = np.zeros(mat["skin"].shape, np.uint8)
    for value, tissue in enumerate(["skin", "bone", "csf", "gray", "white"], start=1):
        m[mat[tissue]] = value
    f = np.zeros(m.shape, np.uint8)
    f[fs_class["ventricles/CSF"]] = 3
    f[fs_class["cortical GM"] | fs_class["subcortical GM"] | fs_class["cerebellar GM"] | fs_class["brainstem"]] = 4
    f[fs_class["white matter"]] = 5
    return m, f


def figure(reference, mat, fs_class, brain_mat, brain_fs, out_png, subject):
    data = reference.get_fdata()
    m, f = overview_label_maps(mat, fs_class)
    agree = np.zeros(m.shape, np.uint8)
    agree[brain_mat & brain_fs] = 1
    agree[brain_mat & ~brain_fs] = 2
    agree[~brain_mat & brain_fs] = 3

    cmap = ListedColormap(PALETTE)
    cmap_agree = ListedColormap([(0, 0, 0, 0), (.4, .8, .4), (.9, .2, .2), (.2, .4, 1)])
    centre = [n // 2 for n in data.shape]
    fig, axes = plt.subplots(2, 3, figsize=(18, 11))
    for row, (axis, idx, plane) in enumerate([(1, centre[1] + 10, "coronal"), (2, centre[2], "axial")]):
        bg = np.take(data, idx, axis=axis).T
        for col, (vol, cm, vmax, title) in enumerate([
            (m, cmap, 5, "MATLAB (SPM + Andy's tools)"),
            (f, cmap, 5, "FreeSurfer (aseg)"),
            (agree, cmap_agree, 3, "brain mask agreement"),
        ]):
            ax = axes[row, col]
            ax.imshow(bg, cmap="gray", origin="lower", vmax=np.percentile(data, 99.5))
            ax.imshow(np.take(vol, idx, axis=axis).T, cmap=cm, origin="lower", vmin=0, vmax=vmax,
                      alpha=.55, interpolation="nearest")
            ax.set_title(f"{plane} {idx}: {title}", fontsize=9)
            ax.set_axis_off()
    axes[0, 0].legend(handles=[Patch(color=c, label=l) for c, l in zip(PALETTE[1:], PALETTE_NAMES)],
                      fontsize=7, loc="lower left")
    axes[0, 2].legend(handles=[Patch(color=c, label=l) for c, l in
                               [((.4, .8, .4), "both"), ((.9, .2, .2), "MATLAB only"), ((.2, .4, 1), "FreeSurfer only")]],
                      fontsize=7, loc="lower left")
    fig.suptitle(f"{subject}: tissue masks, MATLAB pipeline vs FreeSurfer")
    plt.tight_layout()
    plt.savefig(out_png, dpi=85)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--reference", required=True, help="staged NIfTI; defines the comparison grid")
    p.add_argument("--masks", required=True, help="directory with the MATLAB mask_<tissue>.nii files")
    p.add_argument("--aseg", required=True, help="FreeSurfer aseg.mgz")
    p.add_argument("--head", help="FreeSurfer head mask (seghead.mgz)")
    p.add_argument("--out-png", required=True)
    p.add_argument("--out-json", required=True)
    args = p.parse_args()

    reference = nib.load(args.reference)
    voxel_ml = abs(np.linalg.det(reference.affine[:3, :3])) / 1000.0
    subject = os.path.basename(args.reference).split(".")[0]

    mat = {t: nib.load(os.path.join(args.masks, f"mask_{t}.nii")).get_fdata() > 0 for t in MATLAB_TISSUES}
    aseg = resample(args.aseg, reference)
    fs_class = {name: np.isin(aseg, labels) for name, labels in ASEG.items()}
    known = np.concatenate(list(ASEG.values()))
    fs_class["other brain labels"] = (aseg > 0) & ~np.isin(aseg, known)

    brain_mat = mat["gray"] | mat["white"]
    brain_fs = aseg > 0
    head_mat = brain_mat | mat["csf"] | mat["bone"] | mat["skin"]
    head_fs = resample(args.head, reference) > 0 if args.head else None

    # comparable classes: (name, MATLAB mask, FreeSurfer mask)
    pairs = [
        ("white matter", mat["white"], fs_class["white matter"]),
        ("grey matter (all)", mat["gray"],
         fs_class["cortical GM"] | fs_class["subcortical GM"] | fs_class["cerebellar GM"] | fs_class["brainstem"]),
        ("grey matter (cortex only)", mat["gray"], fs_class["cortical GM"]),
        ("CSF", mat["csf"], fs_class["ventricles/CSF"]),
        ("brain (GM+WM)", brain_mat, brain_fs),
    ]
    if head_fs is not None:
        pairs.append(("head", head_mat, head_fs))

    report = {
        "subject": subject,
        "voxel_volume_mm3": round(voxel_ml * 1000, 4),
        "note": "FreeSurfer has no skull, scalp or subarachnoid CSF class; its CSF label covers the ventricles only.",
        "matlab_only_classes_ml": {t: round(float(mat[t].sum()) * voxel_ml, 1) for t in ["bone", "skin"]},
        "classes": {},
        "freesurfer_class_as_matlab_percent": {},
    }
    print(f"{'class':26s} {'MATLAB [ml]':>12s} {'FreeSurfer [ml]':>16s} {'Dice':>6s} {'MAT only':>9s} {'FS only':>9s}")
    for name, a, b in pairs:
        entry = {
            "matlab_ml": round(float(a.sum()) * voxel_ml, 1),
            "freesurfer_ml": round(float(b.sum()) * voxel_ml, 1),
            "dice": round(dice(a, b), 4),
            "matlab_only_ml": round(float((a & ~b).sum()) * voxel_ml, 1),
            "freesurfer_only_ml": round(float((b & ~a).sum()) * voxel_ml, 1),
        }
        report["classes"][name] = entry
        print(f"{name:26s} {entry['matlab_ml']:12.0f} {entry['freesurfer_ml']:16.0f} "
              f"{entry['dice']:6.3f} {entry['matlab_only_ml']:9.0f} {entry['freesurfer_only_ml']:9.0f}")

    print(f"\nFreeSurfer has no equivalent of the MATLAB skull ({report['matlab_only_classes_ml']['bone']:.0f} ml) "
          f"and scalp ({report['matlab_only_classes_ml']['skin']:.0f} ml) classes.")
    print(f"\n{'FreeSurfer class':22s} {'ml':>6s}  " + " ".join(f"{t:>6s}" for t in MATLAB_TISSUES) + f" {'none':>6s}")
    for name, m in fs_class.items():
        n = int(m.sum())
        if not n:
            continue
        frac = {t: round(float((m & mat[t]).sum()) / n * 100, 1) for t in MATLAB_TISSUES}
        frac["unlabelled"] = round(100.0 - sum(frac.values()), 1)
        report["freesurfer_class_as_matlab_percent"][name] = {"volume_ml": round(n * voxel_ml, 1), **frac}
        print(f"{name:22s} {n * voxel_ml:6.0f}  " + " ".join(f"{frac[t]:6.1f}" for t in MATLAB_TISSUES)
              + f" {frac['unlabelled']:6.1f}")

    with open(args.out_json, "w") as fh:
        json.dump(report, fh, indent=2)
    figure(reference, mat, fs_class, brain_mat, brain_fs, args.out_png, subject)
    print(f"\nwrote {args.out_png} and {args.out_json}")


if __name__ == "__main__":
    main()
