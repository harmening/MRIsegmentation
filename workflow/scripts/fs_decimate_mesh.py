#!/usr/bin/env python3
"""Decimate the FreeSurfer cortex mesh so that every vertex keeps voxel support.

    fs_decimate_mesh.py --mesh cortex_pial.obj --out-mesh cortex_pial_25000.obj
                        --reference T1.nii --ribbon ribbon.mgz [options]

The high-resolution pial/white surface of recon-all (~320k vertices) is too
dense for image reconstruction, but plain decimation leaves vertices that no
longer have any voxel of the cortical ribbon mapped onto them - those vertices
cannot be reconstructed. With `--ribbon` (and `--reference`) this decimates in
several steps and, after every step, maps the ribbon voxels onto the vertices:
vertices without voxels ("problematic vertices") and their surroundings get a
low quality score, which steers the next decimation step away from them.

Without `--ribbon` the surface is decimated in one step per entry of
`--nvertices`, without the voxel feedback. That is what the scalp needs: it
carries optodes and bounds the head, but nothing is reconstructed on it, so no
volume has to keep its vertices supported.

With `--annot-lh/--annot-rh` (FreeSurfer .annot files whose vertex order matches
the lh-then-rh order of the mesh, as written by fs_parcels.py) the mapping is
parcel-aware: a voxel is only assigned to vertices of its own parcel, so no
signal leaks across parcel borders. The parcel of each vertex is carried through
the decimation and written to `--out-parcels` in the same format as the MATLAB
pipeline's parcels.json.

`--out-vertexinfo` additionally keeps what the last decimation step computed
anyway - the parcel per vertex and the voxel-to-vertex mapping - so that
make_headmodel.py can build cedalion's voxel_to_vertex_brain matrix and the
per-vertex parcel labels without mapping the voxels a second time.

All meshes are in 1-based voxel indices of `--reference` (the convention of
fs_surfaces_to_voxel.py and the MATLAB exports), and the ribbon is resampled to
that grid, so the output can be compared with the other surfaces directly.

Needs cedalion (with pymeshlab), trimesh, nibabel, scipy, xarray, numpy.
"""
import argparse
import json
import os
import re
import sys

# cedalion's voxel mapping prints a tqdm bar per chunk/parcel, which would fill
# the Snakemake log with thousands of lines
os.environ.setdefault("TQDM_DISABLE", "1")

import cedalion  # noqa: E402
import cedalion.dataclasses as cdc  # noqa: E402
import nibabel as nib  # noqa: E402
import numpy as np  # noqa: E402
import trimesh  # noqa: E402
from cedalion.geometry import meshing  # noqa: E402
from cedalion.geometry import segmentation as segm  # noqa: E402
from scipy.spatial import KDTree  # noqa: E402

# the medial wall / background parcels carry no signal and are not required to
# keep voxel support
BACKGROUND_PREFIX = "Background"


def read_obj(path):
    """Vertices (1-based voxel coords) and 0-based faces of an OBJ file."""
    verts, faces = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith("v "):
                verts.append(line.split()[1:4])
            elif line.startswith("f "):
                faces.append([tok.split("/")[0] for tok in line.split()[1:4]])
    return np.array(verts, float), np.array(faces, int) - 1


def write_obj(path, vertices, faces):
    with open(path, "w") as fh:
        for v in vertices:
            fh.write(f"v {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}\n")
        for f in faces:
            fh.write(f"f {f[0] + 1} {f[1] + 1} {f[2] + 1}\n")


def parcel_name(label, hemi):
    """Parcel name in cedalion's spelling: <name>_LH / <name>_RH.

    17Networks_LH_VisCent_ExStr_1 -> VisCent_ExStr_1_LH. Names without a
    hemisphere (the FreeSurfer medial wall) get the suffix of the hemisphere
    they belong to, so that the two medial walls stay separate parcels and
    match the names cedalion skips when mapping voxels.
    """
    if match := re.match("17Networks_([L,R]H)_(.+)", label):
        return f"{match.group(2)}_{match.group(1)}"
    if label.endswith(("_LH", "_RH")):
        return label
    return f"{label}_{hemi}"


def read_parcels(annot_lh, annot_rh, nvertices):
    """Per-vertex parcel names and colours for the lh-then-rh concatenated mesh."""
    lh_labels, lh_colors, lh_names = nib.freesurfer.read_annot(annot_lh, orig_ids=False)
    rh_labels, rh_colors, rh_names = nib.freesurfer.read_annot(annot_rh, orig_ids=False)
    if len(lh_labels) + len(rh_labels) != nvertices:
        sys.exit(
            f"annotation vertex count {len(lh_labels)} + {len(rh_labels)} does not match "
            f"the mesh ({nvertices}); the .annot files must belong to the surfaces the mesh was built from"
        )
    # labels index into the per-hemisphere name tables; offset the right hemisphere
    all_labels = np.concatenate((lh_labels, rh_labels + len(lh_names)))
    all_colors = np.vstack((lh_colors, rh_colors))
    names = {}
    for i, n in enumerate(lh_names):
        names[i] = parcel_name(n.decode("utf8"), "LH")
    for i, n in enumerate(rh_names):
        names[len(lh_names) + i] = parcel_name(n.decode("utf8"), "RH")
    colors = {names[i]: all_colors[i][:3] for i in names}
    # unlabelled vertices (label -1) belong to the background of their hemisphere
    all_labels = np.where(
        all_labels < 0, np.r_[np.zeros(len(lh_labels), int), np.full(len(rh_labels), len(lh_names))], all_labels
    )
    return np.asarray([names[i] for i in all_labels]), colors


def resample_to_reference(volume_file, reference):
    """Nearest-neighbour resampling of an MGZ/NIfTI volume onto the reference grid."""
    vol = nib.load(volume_file)
    data = np.asarray(vol.dataobj)
    # reference voxel -> world -> volume voxel
    m = np.linalg.inv(vol.affine) @ reference.affine
    grid = np.indices(reference.shape[:3]).reshape(3, -1)
    idx = np.round(m[:3, :3] @ grid + m[:3, 3:]).astype(int)
    inside = np.all((idx >= 0) & (idx < np.array(data.shape[:3])[:, None]), axis=0)
    out = np.zeros(grid.shape[1], data.dtype)
    out[inside] = data[idx[0, inside], idx[1, inside], idx[2, inside]]
    return out.reshape(reference.shape[:3])


def cortical_cell_coords(ribbon_file, reference, mesh_vertices, vertex_parcels):
    """Centres of the cortical ribbon voxels, in 1-based voxel coords of the reference.

    When the surface carries parcels, each voxel inherits the parcel of its
    nearest high-resolution vertex, which is what makes the parcel-aware voxel
    to vertex mapping possible without a separate volumetric parcellation.

    Returns the coordinates, the voxels' flat indices into the full reference
    grid (C order, as np.flatnonzero and cedalion's own voxel-to-vertex map use
    them) and that grid's shape, so that a mapping onto this voxel subset can be
    expressed again over the whole volume.
    """
    ribbon = resample_to_reference(ribbon_file, reference) > 0
    print(f"cortical ribbon: {ribbon.sum()} voxels on the {ribbon.shape} reference grid")

    flat_index = np.flatnonzero(ribbon.ravel())
    coords = segm.cell_coordinates(ribbon, flat=True)
    coords = coords.isel(label=flat_index)
    # cedalion indexes the voxel list positionally, so the label coordinate has
    # to stay 0..N-1 after dropping the voxels outside the ribbon
    coords = coords.assign_coords(label=np.arange(len(coords)))
    # cell_coordinates returns 0-based indices; the meshes use 1-based voxel centres
    coords = coords + 1.0

    if vertex_parcels is not None:
        _, nearest = KDTree(mesh_vertices).query(coords.pint.dequantify().values, workers=-1)
        coords = coords.assign_coords({"parcel": ("label", vertex_parcels[nearest])})
    return coords, flat_index, ribbon.shape


def map_voxels(surface, cell_coords, voxel_stealing, has_parcels):
    """Map the voxels onto the surface and report the vertices without support.

    Returns the vertex index per voxel (-1 where a voxel is not mapped) and a
    mask of the vertices that no voxel is mapped to, background parcels
    excluded - those carry no signal and need no voxel.
    """
    if has_parcels:
        v2v, voxel_count = meshing.parcel_aware_voxels_to_vertices_map(
            surface, cell_coords, voxel_stealing=voxel_stealing
        )
        is_background = surface.vertices.parcel.str.startswith(BACKGROUND_PREFIX).values
    else:
        v2v, voxel_count = meshing.map_voxels_to_vertices(
            surface, cell_coords.pint.dequantify().values
        )
        is_background = np.zeros(surface.nvertices, bool)
    bad = (np.asarray(voxel_count) == 0) & ~is_background
    return np.asarray(v2v), bad


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mesh", required=True, help="high resolution cortex surface (OBJ, 1-based voxel coords)")
    p.add_argument("--reference", help="NIfTI defining the voxel grid of the mesh (with --ribbon)")
    p.add_argument("--ribbon", help="FreeSurfer ribbon.mgz; enables the voxel-support feedback loop")
    p.add_argument("--annot-lh", help="left hemisphere .annot matching the mesh vertex order")
    p.add_argument("--annot-rh", help="right hemisphere .annot matching the mesh vertex order")
    p.add_argument("--out-mesh", required=True)
    p.add_argument("--out-parcels", help="parcel metadata of the decimated mesh (JSON)")
    p.add_argument("--out-vertexinfo",
                   help="per-vertex parcels and the final voxel-to-vertex mapping (NPZ), "
                        "for the head model export; needs --ribbon and the annotations")
    p.add_argument("--nvertices", type=int, nargs="+", default=[30000, 28000, 27000, 26000, 25000],
                   help="decimation schedule; the last value is the target vertex count")
    p.add_argument("--voxel-stealing-from", type=int, default=26000,
                   help="enable voxel stealing at and below this vertex count (parcel-aware mapping only)")
    p.add_argument("--quality-radius", type=float, default=15.0,
                   help="distance (voxels) over which a problematic vertex lowers the decimation quality")
    args = p.parse_args()

    if bool(args.ribbon) != bool(args.reference):
        sys.exit("--ribbon and --reference belong together")

    vertices, faces = read_obj(args.mesh)
    print(f"{args.mesh}: {len(vertices)} vertices, {len(faces)} faces")

    vertex_parcels = None
    colors = {}
    if args.annot_lh and args.annot_rh:
        vertex_parcels, colors = read_parcels(args.annot_lh, args.annot_rh, len(vertices))
        print(f"parcels: {len(set(vertex_parcels))} distinct labels")
    elif args.annot_lh or args.annot_rh:
        sys.exit("--annot-lh and --annot-rh must be given together")

    cell_coords = voxel_flat_index = grid_shape = None
    if args.ribbon:
        cell_coords, voxel_flat_index, grid_shape = cortical_cell_coords(
            args.ribbon, nib.load(args.reference), vertices, vertex_parcels
        )

    # process=False keeps the vertex order, which the parcel labels refer to
    mesh = trimesh.Trimesh(vertices=vertices, faces=faces, process=False)
    surface = cdc.TrimeshSurface(
        mesh,
        crs="ijk",
        units=cedalion.units(None).units,
        vertex_coords={"parcel": vertex_parcels} if vertex_parcels is not None else {},
    )

    quality = None
    v2v = None
    for target in args.nvertices:
        print(f"\ndecimating {surface.nvertices} -> {target}")
        surface = meshing.decimate_mesh(surface, target, quality)

        if cell_coords is None:
            print(f"  {surface.nvertices} vertices")
            continue

        v2v, bad = map_voxels(
            surface, cell_coords,
            voxel_stealing=target <= args.voxel_stealing_from,
            has_parcels=vertex_parcels is not None,
        )
        print(f"  {surface.nvertices} vertices, {bad.sum()} without voxel support")

        if bad.any():
            # vertices close to an unsupported one keep a low quality, so the next
            # step decimates elsewhere
            coords = surface.vertices.pint.dequantify().values
            dists, _ = KDTree(coords[bad]).query(coords, workers=-1)
            quality = np.clip(dists / args.quality_radius, 0.0, 1.0)
        else:
            quality = None

    out_vertices = surface.vertices.pint.dequantify().values
    write_obj(args.out_mesh, out_vertices, surface.mesh.faces)
    print(f"\nwrote {args.out_mesh}: {surface.nvertices} vertices, {surface.nfaces} faces")

    if args.out_parcels:
        if vertex_parcels is None:
            sys.exit("--out-parcels needs --annot-lh/--annot-rh")
        parcels = surface.vertices.parcel.values
        entries = []
        for name in sorted(set(parcels)):
            idx = np.flatnonzero(parcels == name)
            entries.append({
                "Vertices": idx.tolist(),  # 0-based, as in the MATLAB pipeline's parcels.json
                "Color": (np.asarray(colors[name], float) / 255.0).tolist(),
                "Label": name,
                "Region": "RH" if name.endswith("_RH") else "LH",
            })
        with open(args.out_parcels, "w") as fh:
            json.dump(entries, fh, indent=1)
        print(f"wrote {args.out_parcels}: {len(entries)} parcels")

    if args.out_vertexinfo:
        if vertex_parcels is None or v2v is None:
            sys.exit("--out-vertexinfo needs --ribbon/--reference and --annot-lh/--annot-rh")
        # the mapping of the last decimation step, i.e. of the surface just
        # written; voxel_to_vertex and voxel_flat_index are aligned, so the head
        # model can express it over the whole voxel grid again
        np.savez_compressed(
            args.out_vertexinfo,
            parcel=surface.vertices.parcel.values,
            voxel_to_vertex=v2v,
            voxel_flat_index=voxel_flat_index,
            grid_shape=np.asarray(grid_shape),
        )
        print(f"wrote {args.out_vertexinfo}: {np.sum(v2v >= 0)} of {len(v2v)} voxels mapped")


if __name__ == "__main__":
    main()
