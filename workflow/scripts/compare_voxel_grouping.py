#!/usr/bin/env python3
"""QC: what changes when the brain voxels come from SPM instead of FreeSurfer.

    compare_voxel_grouping.py --reference T1.nii --surface cortex_pial_25000.obj
        --fs-subject <SUBJECTS_DIR>/<subj> --grey mask_gray.nii --atlas-dir Schaefer2018_Parcellations
        --out-png voxel_grouping.png --out-json voxel_grouping.json

Image reconstruction groups brain voxels onto the vertices of the cortex
surface (cedalion's voxel-to-vertex map), and every vertex carries a parcel.
The surface and its parcels come from FreeSurfer, but the voxels could come
from either pipeline, so this compares two voxel sets on the *same* surface:

    A  FreeSurfer's cortical ribbon (ribbon.mgz, consistent with the surface)
    B  the MATLAB path's grey matter mask (mask_gray.nii)

For a voxel that both contain the result is identical - the parcel follows
from the surface, not from the mask - so the whole difference comes from which
voxels enter. SPM's "grey matter" is not "cortex": it also covers cerebellum,
brainstem, subcortical grey and part of the white matter, none of which belong
to a cortical parcel but all of which the nearest-vertex mapping assigns to
one. The report quantifies that: voxel-set overlap, what the extra voxels are
anatomically, how far they end up from the vertex they are mapped to, and how
much the per-parcel volumes change.

Schaefer parcels are transferred from fsaverage6 onto the subject's surface by
fs_parcels.py, so no precomputed annotation is required here.

Needs cedalion (trimesh, pymeshlab, rtree), nibabel, numpy, scipy, matplotlib.
"""
import argparse
import json
import os

os.environ.setdefault("TQDM_DISABLE", "1")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import cedalion  # noqa: E402
import cedalion.dataclasses as cdc  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import nibabel as nib  # noqa: E402
import numpy as np  # noqa: E402
import trimesh  # noqa: E402
from cedalion.geometry import meshing  # noqa: E402
from fs_parcels import transfer_parcels  # noqa: E402
from scipy.spatial import KDTree  # noqa: E402

CORTEX_LABELS = [3, 42]  # aseg/ribbon cortical grey matter
# aseg label groups used to describe the voxels SPM adds
ASEG_GROUPS = {
    "cerebellum/brainstem": [8, 47, 7, 46, 16],
    "subcortical grey": [10, 11, 12, 13, 17, 18, 26, 28, 49, 50, 51, 52, 53, 54, 58, 60],
    "cerebral white matter": [2, 41],
    "ventricles/CSF": [4, 5, 14, 15, 24, 43, 44, 72],
}
DISTANCE_THRESHOLDS = (3, 5, 10)


def read_obj(path):
    verts, faces = [], []
    for line in open(path):
        if line.startswith("v "):
            verts.append(line.split()[1:4])
        elif line.startswith("f "):
            faces.append([t.split("/")[0] for t in line.split()[1:4]])
    return np.array(verts, float), np.array(faces, int) - 1


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


def map_set(surface, mask, vertices):
    """Map the voxels of `mask` onto the surface; returns indices and distances."""
    ijk = np.argwhere(mask).astype(float) + 1.0  # 1-based voxel centres, as the meshes use
    v2v, count = meshing.map_voxels_to_vertices(surface, ijk)
    v2v = np.asarray(v2v)
    return ijk, v2v, np.asarray(count), np.linalg.norm(vertices[v2v] - ijk, axis=1)


def figure(reference, data, out_png, subject, ribbon, grey, ijk, dist, is_visual):
    volume_dist = np.zeros(data.shape)
    volume_dist[tuple((ijk - 1).astype(int).T)] = dist
    visual = np.zeros(data.shape, bool)
    visual[tuple((ijk[is_visual] - 1).astype(int).T)] = True

    x = data.shape[0] // 2 + 15
    sl = data[x].T
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    for ax in axes:
        ax.imshow(sl, cmap="gray", origin="lower", vmax=np.percentile(data, 99.5))
        ax.set_axis_off()
    im = axes[0].imshow(np.ma.masked_where(volume_dist[x].T == 0, volume_dist[x].T),
                        cmap="inferno", origin="lower", vmin=0, vmax=30, alpha=.8)
    plt.colorbar(im, ax=axes[0], label="distance to assigned vertex [voxels]")
    axes[0].set_title("SPM grey voxels: distance to the cortex vertex they are mapped to", fontsize=9)
    axes[1].contour(ribbon[x].T, levels=[.5], colors="cyan", linewidths=.6)
    axes[1].contour(grey[x].T, levels=[.5], colors="red", linewidths=.6)
    axes[1].set_title("cyan: FreeSurfer cortical ribbon    red: SPM grey matter", fontsize=9)
    axes[2].imshow(np.ma.masked_where(~visual[x].T, visual[x].T), cmap="autumn", origin="lower", alpha=.8)
    axes[2].set_title("SPM grey voxels assigned to a visual parcel", fontsize=9)
    fig.suptitle(f"{subject}: grouping voxels onto the FreeSurfer cortex surface (sagittal x={x})")
    plt.tight_layout()
    plt.savefig(out_png, dpi=85)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--reference", required=True, help="staged NIfTI; defines the voxel grid")
    p.add_argument("--surface", required=True, help="cortex surface (OBJ, 1-based voxel coords)")
    p.add_argument("--fs-subject", required=True, help="FreeSurfer subject directory")
    p.add_argument("--grey", required=True, help="the MATLAB path's mask_gray.nii")
    p.add_argument("--atlas-dir", required=True, help="directory with the fsaverage6 Schaefer assets")
    p.add_argument("--atlas", default="Schaefer17_600")
    p.add_argument("--out-png", required=True)
    p.add_argument("--out-json", required=True)
    args = p.parse_args()

    reference = nib.load(args.reference)
    data = reference.get_fdata()
    voxel_ml = abs(np.linalg.det(reference.affine[:3, :3])) / 1000.0
    subject = os.path.basename(args.reference).split(".")[0]

    vertices, faces = read_obj(args.surface)
    counts = [len(nib.freesurfer.read_geometry(f"{args.fs_subject}/surf/{h}.pial")[0]) for h in ("lh", "rh")]
    # the surface may be decimated; its vertices are a subset of the original ones
    original, _ = read_obj(os.path.join(os.path.dirname(args.surface), "cortex_pial.obj"))
    parcels_original = transfer_parcels(args.fs_subject, args.atlas_dir, args.atlas, counts)
    _, nearest = KDTree(original).query(vertices, workers=-1)
    parcels = parcels_original[nearest]
    is_cortex_vertex = ~np.char.startswith(parcels.astype(str), "Background")

    surface = cdc.TrimeshSurface(
        trimesh.Trimesh(vertices=vertices, faces=faces, process=False),
        crs="ijk", units=cedalion.units(None).units, vertex_coords={"parcel": parcels},
    )

    aseg = resample(f"{args.fs_subject}/mri/aseg.mgz", reference)
    ribbon = np.isin(resample(f"{args.fs_subject}/mri/ribbon.mgz", reference), CORTEX_LABELS)
    grey = nib.load(args.grey).get_fdata() > 0

    report = {"subject": subject, "atlas": args.atlas, "surface": os.path.basename(args.surface),
              "n_vertices": int(len(vertices)), "n_cortex_vertices": int(is_cortex_vertex.sum())}

    overlap = int((ribbon & grey).sum())
    report["voxel_sets"] = {
        "freesurfer_ribbon_ml": round(float(ribbon.sum()) * voxel_ml, 1),
        "spm_grey_ml": round(float(grey.sum()) * voxel_ml, 1),
        "dice": round(2 * overlap / (int(ribbon.sum()) + int(grey.sum())), 4),
        "spm_only_ml": round(float((grey & ~ribbon).sum()) * voxel_ml, 1),
        "ribbon_only_ml": round(float((ribbon & ~grey).sum()) * voxel_ml, 1),
    }
    print(f"voxel sets: FreeSurfer ribbon {report['voxel_sets']['freesurfer_ribbon_ml']:.0f} ml, "
          f"SPM grey {report['voxel_sets']['spm_grey_ml']:.0f} ml, Dice {report['voxel_sets']['dice']:.3f}")

    only = grey & ~ribbon
    groups = {name: np.isin(aseg, labels) for name, labels in ASEG_GROUPS.items()}
    groups["outside the FreeSurfer brain"] = aseg == 0
    report["spm_only_composition_percent"] = {}
    print("\nthe grey-matter voxels SPM adds on top of the ribbon are, in the FreeSurfer segmentation:")
    for name, m in groups.items():
        share = float((only & m).sum()) / max(int(only.sum()), 1) * 100
        report["spm_only_composition_percent"][name] = round(share, 1)
        print(f"   {name:30s} {float((only & m).sum()) * voxel_ml:6.0f} ml ({share:4.1f} %)")

    mapped = {}
    report["mapping"] = {}
    for key, mask in (("freesurfer_ribbon", ribbon), ("spm_grey", grey)):
        ijk, v2v, count, dist = map_set(surface, mask, vertices)
        mapped[key] = (ijk, v2v, dist)
        unsupported = int(((count == 0) & is_cortex_vertex).sum())
        entry = {
            "n_voxels": int(len(ijk)),
            "unsupported_cortex_vertices": unsupported,
            "unsupported_percent": round(unsupported / int(is_cortex_vertex.sum()) * 100, 2),
            "distance_voxels": {k: round(float(v), 2) for k, v in
                                zip(("median", "p95", "p99", "max"),
                                    (np.median(dist), *np.percentile(dist, [95, 99]), dist.max()))},
            "beyond_threshold_percent": {str(t): round(float((dist > t).mean()) * 100, 1)
                                         for t in DISTANCE_THRESHOLDS},
        }
        report["mapping"][key] = entry
        print(f"\n{key}: {entry['n_voxels']} voxels, {unsupported} cortex vertices without voxels "
              f"({entry['unsupported_percent']:.1f} %)")
        print(f"   voxel-to-vertex distance [voxels]: " +
              ", ".join(f"{k} {v}" for k, v in entry["distance_voxels"].items()))
        print("   beyond: " + ", ".join(f">{t}: {entry['beyond_threshold_percent'][str(t)]} %"
                                        for t in DISTANCE_THRESHOLDS))

    # per-parcel volumes under both voxel sets
    names = sorted(set(parcels[is_cortex_vertex]))
    index = {n: i for i, n in enumerate(names)}
    volumes = {}
    for key, (_, v2v, _) in mapped.items():
        v = np.zeros(len(names))
        for n, c in zip(*np.unique(parcels[v2v], return_counts=True)):
            if n in index:
                v[index[n]] = c * voxel_ml
        volumes[key] = v
    a, b = volumes["freesurfer_ribbon"], volumes["spm_grey"]
    ok = a > 0
    rel = (b[ok] - a[ok]) / a[ok] * 100
    order = np.argsort(rel)
    report["parcel_volumes"] = {
        "n_parcels": int(ok.sum()),
        "total_ribbon_ml": round(float(a.sum()), 1),
        "total_spm_ml": round(float(b.sum()), 1),
        "total_change_percent": round(float((b.sum() - a.sum()) / a.sum() * 100), 1),
        "relative_change_percent": {k: round(float(v), 1) for k, v in
                                    zip(("p5", "p25", "median", "p75", "p95"),
                                        np.percentile(rel, [5, 25, 50, 75, 95]))},
        "largest_increase": [{"parcel": str(np.asarray(names)[ok][i]), "ribbon_ml": round(float(a[ok][i]), 2),
                              "spm_ml": round(float(b[ok][i]), 2), "change_percent": round(float(rel[i]), 0)}
                             for i in order[::-1][:10]],
        "largest_decrease": [{"parcel": str(np.asarray(names)[ok][i]), "ribbon_ml": round(float(a[ok][i]), 2),
                              "spm_ml": round(float(b[ok][i]), 2), "change_percent": round(float(rel[i]), 0)}
                             for i in order[:10]],
    }
    pv = report["parcel_volumes"]
    print(f"\nper-parcel volumes ({pv['n_parcels']} cortical parcels), SPM grey vs FreeSurfer ribbon:")
    print(f"   total {pv['total_spm_ml']:.0f} vs {pv['total_ribbon_ml']:.0f} ml ({pv['total_change_percent']:+.0f} %)")
    print(f"   relative change: median {pv['relative_change_percent']['median']:+.0f} %, "
          f"IQR {pv['relative_change_percent']['p25']:+.0f}..{pv['relative_change_percent']['p75']:+.0f} %, "
          f"5-95 % {pv['relative_change_percent']['p5']:+.0f}..{pv['relative_change_percent']['p95']:+.0f} %")
    for e in pv["largest_increase"][:3]:
        print(f"     {e['parcel']:34s} {e['ribbon_ml']:6.1f} -> {e['spm_ml']:6.1f} ml ({e['change_percent']:+.0f} %)")

    # how much of the "visual cortex" signal would actually be cerebellum
    ijk, v2v, dist = mapped["spm_grey"]
    is_visual = np.char.startswith(parcels[v2v].astype(str), "Vis")
    visual = np.zeros(data.shape, bool)
    visual[tuple((ijk[is_visual] - 1).astype(int).T)] = True
    cereb = np.isin(aseg, ASEG_GROUPS["cerebellum/brainstem"])
    report["visual_parcels"] = {
        "spm_voxels_ml": round(float(visual.sum()) * voxel_ml, 1),
        "of_which_cerebellum_brainstem_ml": round(float((visual & cereb).sum()) * voxel_ml, 1),
    }
    print(f"\nSPM voxels assigned to visual parcels: {report['visual_parcels']['spm_voxels_ml']:.0f} ml, "
          f"of which {report['visual_parcels']['of_which_cerebellum_brainstem_ml']:.0f} ml are "
          f"cerebellum/brainstem")

    with open(args.out_json, "w") as fh:
        json.dump(report, fh, indent=2)
    figure(reference, data, args.out_png, subject, ribbon, grey, ijk, dist, is_visual)
    print(f"\nwrote {args.out_png} and {args.out_json}")


if __name__ == "__main__":
    main()
