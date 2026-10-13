#!/usr/bin/env python3
"""QC: the constructed 10-10 landmarks against a reference set.

    compare_landmarks.py --landmarks headmodel/landmarks.mrk.json
        --scalp headmodel/mask_scalp.obj --reference T1.nii
        --out-png landmarks.png --out-json landmarks.json
        [--reference-landmarks fieldtrip|<file.mrk.json>]

Two things are checked. First, an invariant: the landmarks are constructed by
walking the scalp surface, so every one of them has to lie on it - a nonzero
distance would mean the builder and the shipped surface disagree. Second, with
`--reference-landmarks`, the deviation from an independently defined set, which
is what tells you whether the four input fiducials were sensible.

`fieldtrip` fetches FieldTrip's standard_1005 (defined on the Colin27 head)
through cedalion, i.e. it is downloaded and cached on first use rather than
shipped with this workflow. Anything else is read as a markup file. Without a
reference only the invariant is reported, and the deviation panels are left out.

Needs cedalion, nibabel, numpy, matplotlib, trimesh.
"""
import argparse
import json
import os

os.environ.setdefault("TQDM_DISABLE", "1")

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import nibabel as nib  # noqa: E402
import numpy as np  # noqa: E402
import trimesh  # noqa: E402

FIDUCIALS = ["Nz", "Iz", "LPA", "RPA", "Cz"]
VIEWS = [("sagittal", 1, 2), ("coronal", 0, 2), ("axial", 0, 1)]
AXIS_NAMES = ["x (R)", "y (A)", "z (S)"]


def read_markups(path):
    """Label -> position of a 3D Slicer markup file."""
    with open(path) as fh:
        markups = json.load(fh)["markups"]
    out = {}
    for m in markups:
        for cp in m["controlPoints"]:
            orientation = np.asarray(cp["orientation"], float).reshape(3, 3)
            out[cp["label"]] = np.asarray(cp["position"], float) @ orientation
    return out


def reference_landmarks(spec):
    """The reference set, fetched through cedalion or read from a markup file."""
    if spec == "fieldtrip":
        # downloaded and cached by cedalion on first use, not shipped here
        import cedalion.data

        lm = cedalion.data.get_fieldtrip_colin27_landmarks()
        pos = lm.pint.to("mm").pint.dequantify().values
        return {str(x): np.asarray(p) for x, p in zip(lm.label.values, pos)}, "FieldTrip standard_1005"
    return read_markups(spec), os.path.basename(spec)


def figure(out_png, subject, scalp_ras, built, ref, ref_name, on_scalp, deviations):
    # three views plus a histogram, which is the deviation from the reference
    # when there is one and the on-surface check otherwise
    fig, axes = plt.subplots(1, 4, figsize=(24, 6.0))
    rng = np.random.default_rng(0)
    thin = scalp_ras[rng.choice(len(scalp_ras), size=min(4000, len(scalp_ras)), replace=False)]

    labels = sorted(built)
    for ax, (name, a, b) in zip(axes[:3], VIEWS):
        ax.scatter(thin[:, a], thin[:, b], s=1, c="0.85", linewidths=0)
        if ref is not None:
            common = [k for k in labels if k in ref]
            r = np.array([ref[k] for k in common])
            ax.scatter(r[:, a], r[:, b], s=26, facecolors="none", edgecolors="goldenrod",
                       label=f"reference ({ref_name})")
        p = np.array([built[k] for k in labels])
        ax.scatter(p[:, a], p[:, b], s=12, c="tab:green", label="constructed")
        f = np.array([built[k] for k in FIDUCIALS if k in built])
        ax.scatter(f[:, a], f[:, b], s=40, c="tab:red", label="fiducials")
        for k in FIDUCIALS:
            if k in built:
                ax.annotate(k, (built[k][a], built[k][b]), fontsize=7,
                            textcoords="offset points", xytext=(4, 3))
        ax.set_title(f"{subject} - {name}")
        ax.set_xlabel(AXIS_NAMES[a])
        ax.set_ylabel(AXIS_NAMES[b])
        ax.set_aspect("equal")
    axes[0].legend(loc="lower left", fontsize=7)

    ax = axes[-1]
    if ref is None:
        ax.hist(on_scalp, bins=np.arange(0, max(1.0, on_scalp.max()) + 0.05, 0.05), color="tab:green")
        ax.set_xlabel("distance of a landmark from the scalp / mm")
        ax.set_title("on-surface check")
    else:
        d = np.array(list(deviations.values()))
        ax.hist(d, bins=np.arange(0, max(10.5, d.max() + 0.5), 0.5), color="tab:green")
        ax.axvline(np.median(d), c="r", ls=":", label=f"median {np.median(d):.1f} mm")
        ax.set_xlabel(r"$|\vec{x}_{ours} - \vec{x}_{ref}|$ / mm")
        ax.set_title(f"deviation from {ref_name}")
        ax.legend(fontsize=8)
    ax.set_ylabel("# landmarks")

    plt.tight_layout()
    plt.savefig(out_png, dpi=90)
    print(f"wrote {out_png}")


def main():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    p.add_argument("--landmarks", required=True, help="constructed landmarks (.mrk.json)")
    p.add_argument("--scalp", required=True, help="scalp surface (OBJ, 1-based voxel indices)")
    p.add_argument("--reference", required=True, help="NIfTI defining the voxel grid of the mesh")
    p.add_argument("--reference-landmarks", default="",
                   help="'fieldtrip' or a markup file; empty = no comparison")
    p.add_argument("--out-png", required=True)
    p.add_argument("--out-json", required=True)
    p.add_argument("--subject", default="")
    args = p.parse_args()

    built = read_markups(args.landmarks)
    voxel_mesh = trimesh.load(args.scalp, process=False)
    # the scalp is in 1-based voxel indices, the affine maps 0-based voxel centres
    mesh = trimesh.Trimesh(
        vertices=nib.affines.apply_affine(
            np.asarray(nib.load(args.reference).affine, float),
            np.asarray(voxel_mesh.vertices) - 1.0,
        ),
        faces=voxel_mesh.faces,
        process=False,
    )
    scalp_ras = np.asarray(mesh.vertices)
    print(f"{args.landmarks}: {len(built)} landmarks; scalp: {len(scalp_ras)} vertices")

    # invariant: the landmarks were walked along this surface, so they lie on it.
    # Measured to the nearest point of a face, not to the nearest vertex - the
    # builder places them along face intersections, so a vertex distance would
    # only report the mesh's edge spacing.
    labels = sorted(built)
    _, on_scalp, _ = trimesh.proximity.closest_point(
        mesh, np.array([built[k] for k in labels])
    )
    edges = mesh.edges_unique_length
    print(f"distance from the scalp: median {np.median(on_scalp):.3f} mm, "
          f"max {on_scalp.max():.3f} mm (edge length: median {np.median(edges):.2f} mm)")

    report = {
        "subject": args.subject,
        "n_landmarks": len(built),
        "distance_from_scalp_mm": {
            "median": float(np.median(on_scalp)), "max": float(on_scalp.max()),
        },
        "scalp_edge_length_mm": {
            "median": float(np.median(edges)), "max": float(edges.max()),
        },
        "reference": None,
    }

    ref = ref_name = None
    deviations = {}
    if args.reference_landmarks:
        ref, ref_name = reference_landmarks(args.reference_landmarks)
        common = sorted(set(built) & set(ref))
        deviations = {k: float(np.linalg.norm(built[k] - ref[k])) for k in common}
        d = np.array(list(deviations.values()))
        print(f"vs {ref_name}: {len(common)} common labels, median {np.median(d):.2f} mm, "
              f"mean {d.mean():.2f}, max {d.max():.2f} ({common[int(np.argmax(d))]})")
        report["reference"] = {
            "name": ref_name,
            "spec": args.reference_landmarks,
            "n_common": len(common),
            "median_mm": float(np.median(d)),
            "mean_mm": float(d.mean()),
            "max_mm": float(d.max()),
            "max_label": common[int(np.argmax(d))],
            "only_in_ours": sorted(set(built) - set(ref)),
            "only_in_reference": sorted(set(ref) - set(built)),
            "deviation_mm": deviations,
        }
    else:
        print("no reference set given, reporting the on-surface check only")

    figure(args.out_png, args.subject or "subject", scalp_ras, built, ref, ref_name,
           on_scalp, deviations)
    with open(args.out_json, "w") as fh:
        json.dump(report, fh, indent=2)
    print(f"wrote {args.out_json}")


if __name__ == "__main__":
    main()
