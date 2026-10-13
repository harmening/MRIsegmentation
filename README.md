# MRIsegmentation
[![DOI](https://zenodo.org/badge/457134654.svg)](https://zenodo.org/badge/latestdoi/457134654)
### Automatic MRI segmentation pipeline for consistent [FEM](https://en.wikipedia.org/wiki/Finite_element_method) and [BEM](https://en.wikipedia.org/wiki/Boundary_element_method) mesh creation from MRI scans of human heads
**0. <ins>[Install prerequisites and start the pipeline](#prerequisites)</ins>**<br>

**1. <ins>[Translation to ACPC-coordinate system](#translation-to-acpc)<ins>**<br>
**2. <ins>[Tissue segmentation (Scalp, Skull, CSF, White Matter, Gray Matter)](#tissue-segmentation)</ins>**<br>
**3. <ins>[Surface and Volume Mesh construction](#mesh-construction)</ins>**<br>
**4. <ins>[Correction and smoothing of segmentation mistakes](#correction-and-smoothing)</ins>**<br>
**5. <ins>[Alignement of electrodes (and sourcemodels)](#electrode-alignement)</ins>**<br>
**6. <ins>[Cortical sourcemodel creation](#sourcemodel-creation)</ins>**<br>
**7. <ins>[Artefact sourcemodel creation](#hartmuts-artefact-sourcemodel)</ins>**<br>
**8. <ins>[Transformation to fiducial based coordinate systems](#transformation-to-fiducial-based-coordinate-systems)</ins>**<br>
**9. <ins>[Representing the cortical surface in Spherical Harmonics (SPHARM)](#representing-the-cortical-surface-in-spherical-harmonics)</ins>**<br>
**10. <ins>[Labeling with parcellation atlas](#parcellation-labels)</ins>**<br>

**A1. <ins>[Troubleshooting](#troubleshooting)</ins>**<br>
**A2. <ins>[Citing](#citing)</ins>**<br>
**A3. <ins>[References](#references)</ins>**<br>
<br>
<br>


## Prerequisites
### This segmentation pipeline connects these awesome open source projects:
* [Automatic Temporal Registration Algorithm (ATRA)](https://www.nitrc.org/projects/art) [[1]](#ref1)
* [Statistical Parametric Mapping (SPM12)](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) [[2]](#ref2)
* [Andy's tools](https://www.parralab.org/segment/) [[3]](#ref3)
* [iso2mesh](http://iso2mesh.sourceforge.net)	[[4]](#ref4)
* [Computational Anatomy Toolbox (CAT12)](http://www.neuro.uni-jena.de/cat/)	[[5]](#ref5)
* [fieldtrip](http://www.fieldtriptoolbox.org/) [[6]](#ref6)
* [HArtMuT](http://www.hartmut.berlin/) [[7]](#ref7)
* [weighted SPHARM](https://pages.stat.wisc.edu/~mchung/softwares/weighted-SPHARM/) [[8]](#ref8)
* [Open Access Series of Imaging Studies (OASIS)](http://oasis-brains.org/) [[9]](#ref9) <br>
* [Schaefer2018 Brain Parcelalation Atlas](https://github.com/ThomasYeoLab/CBIG/tree/v0.14.3-Update_Yeo2011_Schaefer2018_labelname/stable_projects/brain_parcellation/Schaefer2018_LocalGlobal) [[10]](#ref10) <br>
<br>

<!--- 
[extended Tissue Probability Map (eTPM.nii)](http://neuralengr.com/segment/) [[3]](#ref3)
freesurfer (conversion to nifti)
--->

### How to get up and running:
- Install Matlab and its [NIfTI and ANALYZE tools](https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image).
<!--- 
- Install [fieldtrip](http://www.fieldtriptoolbox.org/download/). [SPM12](https://www.fil.ion.ucl.ac.uk/spm/software/spm12/) should come with it; if not, install it, too.
- Download ATRA for [Linux](https://www.nitrc.org/frs/download.php/10393/atra1.0_LinuxCentOS6.7.tar.gz) / [Mac](https://www.nitrc.org/frs/download.php/10394/atra1.0_macOS.10.12.6.tar.gz) into `/usr/local/art/`.
--->
- Run the following commands in a terminal:<br>

<!--- 
- Install [CAT12](http://www.neuro.uni-jena.de/cat/index.html#DOWNLOAD), unzip and move the cat12 folder into `fieldtrip/external/spm12/toolbox`.
https://neuroimage.usc.edu/brainstorm/Tutorials/SegCAT12
--->
```bash
# Clone repository
$ git clone https://github.com/harmening/MRIsegmentation.git
$ cd MRIsegmentation
# Install fieldtrip
$ git clone https://github.com/fieldtrip/fieldtrip.git
# Install SPM12 
$ wget -P ./ https://www.fil.ion.ucl.ac.uk/spm/download/restricted/eldorado/spm12.zip
$ rm -r ./fieldtrip/external/spm12
$ unzip spm12.zip -d ./fieldtrip/external/
$ cd ./fieldtrip/external/spm12/src
$ make distclean
$ make && make install
$ make external-distclean
$ make external && make external-install
$ cd ../../../../ && rm spm12.zip
# Install CAT12
$ wget -P ./ http://www.neuro.uni-jena.de/cat12/cat12_r2159.zip
$ rm -r ./fieldtrip/external/spm12/toolbox/cat12
$ unzip cat12_r2159.zip -d ./fieldtrip/external/spm12/toolbox/
$ rm cat12_r2159.zip
# Install iso2mesh
$ git clone https://github.com/fangq/iso2mesh.git
# Download Andy's tools
$ wget -P ./ https://www.parralab.org/segment/Huang_et_al_2013.zip
$ unzip Huang_et_al_2013.zip -d ./Huang_et_al_2013
$ rm Huang_et_al_2013.zip
# Minor changes to make new_segment storing the nonlinear warp
$ sed -i.bak 's/warp\.write = \[0 0\]/warp\.write = \[1 1\]/' Huang_et_al_2013/start_seg.m  && rm Huang_et_al_2013/start_seg.m.bak
# Load HArtMuT Individual Warp
$ git clone https://github.com/harmening/hartmut.git HArtMuT
```
- The coordinate systems are defined by the MRI scanner's coordinate systems. For optional [translation into the ACPC-coordinate system](#translation-to-acpc) before segmentation, download [atra for linux](https://www.nitrc.org/frs/download.php/10393/atra1.0_LinuxCentOS6.7.tar.gz/?i_agree=1&release_id=3672) (register at nitrc.org, agree to the license terms), move the downloaded `atra1.0_LinuxCentOS6.7.tar.gz` into `./art` and run the following commands in a terminal:<br>
```bash
# Install atra
$ export ARTHOME="$(pwd)/art"; export PATH=$ARTHOME/bin:$PATH
$ cd $ARTHOME
$ gunzip atra1.0_*.tar.gz
$ tar -xvf atra1.0_*.tar
$ rm atra1.0_*.tar; rm atra1.0_*.tar.gz;
```
- If desired, download the weighted-SPHARM ([`SPHARMsmooth2.m`](https://pages.stat.wisc.edu/~mchung/research/amygdala/SPHARMsmooth2.m), [`SPHARMvectorize.m`](https://pages.stat.wisc.edu/~mchung/softwares/weighted-SPHARM/SPHARMvectorize.m)) into `./weighted-SPHARM` by following the [description of Chung et al.](https://pages.stat.wisc.edu/~mchung/softwares/weighted-SPHARM/).


### How to start segmentation:
- The MRI scan needs to be in NIfTI file format (.nii or .nii.gz).
```bash
# Change parameters in start_segmentation.m (EEG pipeline) and run 
$ matlab -r "start_segmentation(./example_head/T1.nii); exit;"
# For fNIRS/DOT run
$ matlab -r "start_segmentation_cedalion(./example_head/T1.nii); exit;"
```

### Container + Snakemake workflow
Instead of installing everything by hand, build the container image (MATLAB via
[mpm](https://github.com/mathworks-ref-arch/matlab-dockerfile), SPM12, CAT12,
FieldTrip, iso2mesh, Andy's tools, NIfTI tools, HArtMuT, weighted-SPHARM) and
run the pipeline step by step with [Snakemake](https://snakemake.readthedocs.io):
```bash
# 1. Dependency archives (SPM12, CAT12 r2159, Andy's tools) next to the Dockerfile.
#    Optional: put atra1.0_LinuxCentOS6.7.tar.gz there too for ACPC preprocessing.
$ bash docker/fetch-deps.sh
# 2. Build the image (podman or docker; ~12 GB, MATLAB is downloaded by mpm) and
#    convert it for apptainer (~4 GB mriseg.sif, see licensing below).
#    --build-arg MATLAB_RELEASE=R2024b / MATLAB_TOOLBOXES="..." to change defaults.
$ podman build -t mrisegmentation .
$ bash docker/build-sif.sh
# 3. Describe your subjects (tab separated: subject, t1, t2) and settings
$ $EDITOR config/samples.tsv config/config.yaml   # pipeline: cedalion | eeg
# 4. Run
$ pip install snakemake
$ snakemake -n              # dry run: shows the per-subject jobs
$ snakemake --cores 2       # two subjects in parallel
```
Outputs are written to `results/<subject>/` (the same files the entry scripts
produce, next to the staged `<subject>.nii`), MATLAB logs to
`results/logs/<subject>/<step>.log`. Each pipeline function is one rule
(`spm_segment` → `tissue_masks` → `surface_meshes` / `volume_meshes` → `cat12` →
`parcellation` / `mne_export`, plus `align_electrodes`, `hartmut_*`, `sourcemodel`,
`transform_to_ctf`, `spharm` for the EEG pipeline), so a failed or changed step
only reruns from there: `snakemake --cores 1 results/sub01/parcels.json`.

**Tissue probability map.** SPM segments with the extended eTPM of Andy's
tools by default (skull/scalp/air priors down to the neck). On some template
heads it collapses — Colin27 ends up with empty GM/WM/CSF classes and the
workflow now stops at `tissue_masks` with a clear error instead of failing
later in the meshing. For such inputs set `spm_tpm: TPM` (SPM12's own
`tpm/TPM.nii`, same six classes, coarser skull/scalp priors and no neck) in
`config/config.yaml`, or per subject via a `tpm` column in the sample sheet.

**FreeSurfer path (optional).** With `freesurfer.enabled: true` the workflow
also runs `recon-all` (+ `mkheadsurf` for the scalp) per subject in a separate
FreeSurfer container, in parallel to the MATLAB steps, and exports the surfaces
into the same coordinate system as the MATLAB results for direct comparison:
```bash
$ apptainer pull freesurfer.sif docker://freesurfer/freesurfer:7.4.1   # ~9 GB
$ $EDITOR config/config.yaml    # freesurfer.license: path to your license.txt
```
`results/<subject>/fs/<subject>/` is a regular FreeSurfer subject directory
(`SUBJECTS_DIR=results/<subject>/fs`), `results/<subject>/fs/*.obj` are the
scalp (`scalp`, `scalp_smoothed`), `lh/rh.white`, `lh/rh.pial` and combined
`cortex_white` / `cortex_pial` surfaces with vertices as 1-based voxel indices
of the staged NIfTI — the convention of `mask_scalp.obj` / `mask_brain.obj`
(see `workflow/scripts/fs_surfaces_to_voxel.py` for the tkr-RAS → voxel chain).
A T2 in the sample sheet is used for pial refinement (`-T2pial`); inputs taller
than 256 mm get `-cw256` automatically. `mkheadsurf`'s head mask is cleaned
before tessellation: a small morphological opening (`freesurfer.scalp_open_radius`)
removes the thin sheets that ghosting artefacts next to the head otherwise leave
on the scalp surface, and a population prior (`freesurfer.scalp_margin_mm`,
`freesurfer.scalp_dim_below`) drops voxels that lie above the vertex, outside
the MNI305 population head (via the subject's `talairach.xfm`) and dark in the
normalised T1 — structures that touch the scalp and therefore survive any
opening, such as the bar across the top of the Colin27 template. A head that is
merely larger than the template keeps its bright scalp. `scalp_dim_below` must
stay below the intensity of the outer scalp rim, which fades to ~40-70 in the
normalised T1: a higher value shaves the crown instead of only the artefact.
`mri/seghead.orig.mgz` keeps the uncleaned mask for inspection.

With `freesurfer.decimate.enabled` the high-resolution surfaces are
additionally reduced, to `results/<subject>/fs/<source>_<N>.obj`;
`freesurfer.decimate.targets` maps each `fs_export` surface to its decimation
schedule (cortex ~320k -> 25k vertices, scalp ~200k -> 10k by default).
Surfaces listed under `ribbon_support` - the cortex surfaces, on which images
are reconstructed - must keep cortical-ribbon voxels mapped onto every vertex,
which plain decimation does not guarantee, so
`workflow/scripts/fs_decimate_mesh.py` decimates in several steps and, after
each one, maps the ribbon voxels (`mri/ribbon.mgz`, resampled to the NIfTI grid)
onto the vertices: vertices without voxel support and their neighbourhood get a
low quality score that steers the next step away from them. Given surface
annotations (`freesurfer.decimate.annot`) the mapping is parcel-aware and the
per-vertex parcels are carried through the decimation and written next to the
mesh as `<source>_<N>_parcels.json`, in the same format as the MATLAB pipeline's
`parcels.json`, plus `<source>_<N>_vertexinfo.npz` for the head model export
below. recon-all ships no Schaefer annotation, so rule `fs_annot` writes one:
`workflow/scripts/fs_parcels.py` transfers `freesurfer.parcel_atlas` from
fsaverage6 onto the subject by nearest neighbour on `sphere.reg` (the transfer
`mri_surf2surf` would do) into the subject's `label/` directory, with the parcel
names and colours the atlas ships. The scalp carries optodes and bounds the head but nothing is
reconstructed on it, so it is simply decimated in one step (point-to-surface
error at 10k vertices: median 0.02, max 1.9 voxels). This rule runs on the host
rather than in a container: it needs
[cedalion](https://github.com/ibs-lab/cedalion) (`cedalion.geometry.meshing`)
plus trimesh, pymeshlab and rtree in the environment that runs Snakemake, or in
the one given as `freesurfer.decimate.python`. recon-all takes a few hours per subject
(`freesurfer.threads`, default 8, for `-parallel -openmp`). With both paths
present, `results/qc/<subject>/surfaces_matlab_vs_freesurfer.{png,json}` shows
exact cross-sections of the MATLAB scalp / CAT12 cortex and the FreeSurfer
scalp / pial / white surfaces on the T1 plus nearest-vertex distances, and
`masks_matlab_vs_freesurfer.{png,json}` compares the volumetric segmentations:
volumes and Dice of white matter, grey matter, CSF, brain and head, plus a
cross-tabulation of how the MATLAB masks label each FreeSurfer class. Note that
the two segment different things - recon-all has no skull, scalp or
subarachnoid CSF class, and its CSF label covers the ventricles only - so only
brain classes and the head mask are really comparable.

A third QC rule, `voxel_grouping_matlab_vs_freesurfer.{png,json}`, answers what
it costs to mix the two: image reconstruction groups brain voxels onto the
vertices of the cortex surface, and the rule maps both candidate voxel sets -
FreeSurfer's cortical ribbon and the MATLAB path's `mask_gray.nii` - onto the
same decimated surface, with Schaefer parcels transferred from fsaverage6
(`freesurfer.parcel_atlas`). It reports the overlap of the two voxel sets, what
the voxels SPM adds are anatomically, how far they end up from the vertex they
are mapped to, and how the per-parcel volumes change. The first two rules need
numpy, scipy, nibabel and matplotlib; this one also needs cedalion, like
`fs_decimate`.

**cedalion head model.** With `freesurfer.headmodel.enabled` the workflow
assembles `results/<subject>/headmodel/`, a directory in the layout
[cedalion](https://github.com/ibs-lab/cedalion) expects of a head model
(`TwoSurfaceHeadModel.from_surfaces`, `cedalion.data.get_*_headmodel_files`), as
in its published `hm_colin27` / `hm_icbm152` datasets:

| file | content |
| --- | --- |
| `mask_{csf,gray,white,bone,skin}.nii` | tissue masks of the MATLAB path |
| `cortex_pial_high.obj` | unreduced FreeSurfer cortex surface |
| `mask_brain.obj` | the decimated cortex surface (25k vertices) |
| `cortex_pial_inflated.obj` | the same vertices on the inflated cortex |
| `mask_scalp.obj` | the decimated scalp surface (10k vertices) |
| `landmarks.mrk.json` | the 10-10 system, built on `mask_scalp.obj` |
| `brain_vertex_coordinates.csv` | per vertex: parcel, MNI152 position, and the vertex it came from on `cortex_pial_high` |
| `voxel_to_vertex_brain.mtx.gz` | cortical-ribbon voxel → cortex vertex, as a MatrixMarket matrix of one row per voxel of the whole grid |
| `parcel_colors.json` | parcel name → 8-bit RGB |

The voxel-to-vertex matrix is the parcel-aware mapping of the last decimation
step, so it is built from FreeSurfer's cortical ribbon rather than from SPM's
grey plus white matter, which also covers cerebellum, brainstem and subcortical
grey — see the `voxel_grouping` QC above for what that difference costs. For the
`mni152_*` columns, `freesurfer.headmodel.mni_from` (overridable per subject
with an `mni_from` column in the sample sheet) states what the image's own RAS
already is: `talairach` for an individual scan, which recon-all's
`talairach.xfm` registers to MNI305 before cedalion's fixed MNI305 → MNI152
affine; `mni305` for a template distributed in MNI305, such as Colin27, which
needs only that last step (what the published cedalion head models did); and
`mni152` for the ICBM152 templates, whose RAS is MNI152 already, so that
converting them from MNI305 would only add error. The shipped sample sheet marks
the two templates accordingly.

The meshes keep this pipeline's **1-based** voxel indices, which is what cedalion
expects of a head model's surface files: `TwoSurfaceHeadModel.from_surfaces`
shifts them into its own 0-based `ijk` itself (its `surfaces_are_zero_based`
argument). The `mni152_*` columns are computed after that shift, so they describe
the vertex where cedalion puts it. `cortex_pial_inflated.obj` is not a voxel grid
at all — the inflated cortex is its own coordinate system — and the rows of
`voxel_to_vertex_brain.mtx.gz` are 0-based flat voxel indices either way.

**Landmarks.** Rule `landmarks` builds the 10-10 system on `mask_scalp.obj` with
cedalion's `LandmarksBuilder1010`, which walks the scalp between the fiducials
`Nz`, `Iz`, `LPA` and `RPA` ([Oostenveld and Praamstra
2001](https://doi.org/10.1016/S1388-2457(00)00527-7)). Those four are *inputs*,
not derived: each subject points at its own 3D Slicer markup file with a
`landmarks` column in the sample sheet, in the RAS of its staged NIfTI (see
[config/landmarks/README.md](config/landmarks/README.md); files for the two
templates are included). `Cz` is optional, as the builder re-estimates the
cranial vertex anyway. The fiducials are snapped onto the scalp first
(`freesurfer.landmarks.snap`) — literature fiducials are defined on another
head's scalp, so a few mm of offset is normal — and a fiducial further than
`max_snap_mm` aborts the rule, which is what a markup file in the wrong
coordinate system looks like. A file exported from 3D Slicer with its default
`LPS` is converted; an unrecognised `coordinateSystem` is rejected rather than
guessed at.

`results/qc/<subject>/landmarks.{png,json}` checks the result two ways: that
every landmark really lies on the scalp (it is constructed by walking it, so the
distance must be ~0 — measured to the nearest point of a face, not the nearest
vertex), and, where a reference set is configured with a `landmarks_reference`
column, the deviation from it. `fieldtrip` fetches FieldTrip's `standard_1005`
(defined on Colin27, 346 points) through cedalion, i.e. it is downloaded and
cached on first use rather than shipped here; any other value is read as a
markup file, and an empty one skips the comparison.

Expect the landmarks to differ from the ones in cedalion's published head
models by a few mm: those were constructed on cedalion's SPM-derived scalp,
these on the FreeSurfer scalp this directory ships, which is the surface optodes
are later placed on.

**MATLAB licensing inside the container.** By default the image is built in
online-licensing mode (`--build-arg MATLAB_LICMODE=onlinelicensing`), i.e. MATLAB
licenses through your MathWorks account exactly like a desktop install where you
chose "sign in". No credentials go into the image: the workflow runs the `.sif`
with apptainer, which executes the container as *you* with your home directory,
hostname and network, so MATLAB inside finds the sign-in your host MATLAB cached
in `~/.matlab` / `~/.MathWorks` (sign in once on the host with the same release;
Individual licenses work this way without a token). podman/docker do not share
that identity — the cached credentials fail their integrity check there — so for
those build with `--build-arg MATLAB_LICMODE=file` and provide, in your
environment, `MLM_LICENSE_FILE=27000@licenseserver` (network license),
`container.license_file: /path/to/license.lic` (license file) or
`MLM_LICENSE_TOKEN` ([batch licensing token](https://www.mathworks.com/support/batch-tokens.html)).

Sharing your home that way also shares `~/.MathWorks/ServiceHost/<hostname>`, which
besides the cached sign-in holds the runtime state of the MathWorks Service Host — a
daemon that detaches and outlives the MATLAB that started it. Shared, every job and
every MATLAB on the host (e.g. the one the VS Code extension starts) compete for one
daemon whose socket moreover lives in a per-job `TMPDIR` that is gone by the time the
next job connects, which surfaces as `Unable to communicate with required MathWorks
services` (error 5001/5201). Each job therefore gets its own copy of that directory
(~25 kB without the logs) bound over the shared one, and runs in its own PID namespace
so a detached daemon cannot keep the container alive. Nothing is written back, so
refresh your sign-in on the host when it expires.

Unrelated but easy to confuse with a licensing problem: do not run MATLAB under
`ulimit -v`. It caps the virtual address space, `pthread_create` then fails with
`EAGAIN`, and the Service Host aborts in its gRPC server — reported, again, as error
5001/5201.

**Apptainer on Ubuntu ≥ 23.10.** Two one-time system tweaks: install
`squashfuse` so the `.sif` is mounted instead of being extracted (3.8 GB into
`/tmp`) for every job, and give Apptainer permission to create user namespaces
(Ubuntu's apptainer package has no AppArmor profile; without it every run fails
with `Failed to create container process: Operation not permitted`):
```bash
$ sudo apt install squashfuse
$ sudo cp docker/apparmor-apptainer /etc/apparmor.d/apptainer && sudo apparmor_parser -r /etc/apparmor.d/apptainer
```

The image is also usable without Snakemake:
`apptainer exec mriseg.sif mriseg "start_segmentation_cedalion('/data/sub01/T1.nii')"`
runs a whole entry script; `mriseg` accepts any MATLAB statement after setting
up the path (`docker/mriseg_init.m`). With podman use
`podman run --rm --user root -e MLM_LICENSE_FILE -v /data:/data mrisegmentation mriseg "..."`
(rootless podman: container root is your user, so the outputs belong to you;
with docker use `--user $(id -u):$(id -g)`).
<!--- 
# Or alternatively in parallel:
$ python parallel.py ./data/MRI_scan1/T1.nii ./data/MRI_scan2/T1.img ./data/MRI_scan3/T1.nii.gz
--->
<br>
<br>
<br>




## Translation to ACPC
 - Anatomical landmark detection, especially anterior and posterior 
   <img align="right" width="200" src="images/t1_mprage_orion.png"> <br>
 commissures (AC/PC).
 - Reorientation from MRI scanner coordinates into the comparable ACPC <br>
   coordinate system (origin at AC, RAS orientation (x = right, y = anterior, <br>
   z = superior), original size i.e. not normalized to a template)
 - For more details refer to [[1]](#ref1).<br>
<br>

## Tissue segmentation
- The posterior probability of each MRI voxel belonging to a specific tissue type
<img align="right" width="200" src="images/segmentedmri.png"> <br>
  (scalp, skull, cerebrospinal fluid (CSF), gray matter (GM), white matter(WM)) <br>
  is estimated by the "Unified Segmentation" framework of [[2]](#ref2) as implemented <br>
  as *new_segment* in SPM12.
- As prior probability distribution template the extended tissue probability map (Huang_et_al_2013/eTPM.nii) developed in [[3]](#ref3) is used.
- The calculated probability distributions are stored (c1, c2, ..) and converted <br>
  into binary masks for each tissue (mask_scalp.nii, mask_skull.nii, ..), followed <br>
  by smoothing using a Gaussian low-pass filter using *Andys tools* [[3]](#ref3).
- The following automatic correction of obvious morphological errors includes assigning unassigned voxels and reclassifying GM voxels, that are directly adjacent to skull voxels, as CSF voxels in order to avoid potential CSF discontinuities. The final segmented MRI is stored as segmentedmri.mat.<br>
<br>

## Mesh construction
- Before meshing, labelled voxels that cannot belong to the head are dropped
  (`restrict_to_head.m`): voxels whose air prior in the tissue probability map,
  sampled in native space through SPM's deformation field `iy_*.nii`, exceeds
  0.99, and everything outside the largest connected non-air region. This
  removes ghosting/averaging artefacts next to the head that SPM labels as
  tissue when the background is very clean (e.g. the lateral smears of the
  ICBM152 template), which otherwise end up as fins on the scalp mesh.
- Boundary surface meshes per tissue are triangulated with *cgalsurf* as 
<img align="right" width="250" src="images/volmesh.png"> <br>
  implemented in the fieldtrip toolbox [[6]](#ref6). From these meshes (bnd5.mat) <br>
  and (bnd4.mat, merged GM and WM to one cortex mesh) corresponding <br>
  5-shell and 4-shell Boundary Element Method (BEM) models can be built.
- A volumetric tetrahedral mesh is constructed using vol2mesh [[4]](#ref4). <br>
  Using *cgalmesh* it creates high-quality Delaunay tetraeder (see <br>
  plotted tissued voxels to the right), that can be used for Finite <br>
  Element Method (FEM) modelling. Nodes, faces, elements and voxel <br>
  labels are stored as mesh5.mat.<br>
<br>

## Correction and smoothing 
- Tests on the [OASIS1 database](http://oasis-brains.org/) of 416 human heads [[9]](#ref9) revealed, that 
<img align="right" width="250" src="images/cortex_correct.gif"> <br>
  10 out of 416 segmentations needed further correction of the cortical <br>
  surface mesh only (see segmentation error and correction to the right).
- The correction is done by recalculating abnormal vertex positions <br>
  (as in example to the right) relative to their neighboring mesh nodes.
- For an extra smooth cortex (GM+WM) surface a further Laplace flow <br>
  mesh smoothing is applied.<br>
<br>

## Electrode alignement
- Nonlinear alignment of [standard electrode positions](http://www.fieldtriptoolbox.org/template/electrode/) by using SPMs nonlinear,
<img align="right" width="250" src="images/elec_aligned.png"> <br>
voxel-wise transformation mapping from the eTPM (in MNI space) to the individual head coordinates.
- Refinement of electrode positions by projecting them directly onto the segmented scalp mesh. This small projection corrects for minimal mesh approximation errors and ensures that all electrodes are perfectly placed on the outermost mesh and thus can be used directly for BEM modelling.
- Accurate nonlinear warping of [fiducials (nasion (NAS), left and right pre-auricular (LPA/RPA))](http://www.fieldtriptoolbox.org/template/) gained from manual expert selection on a high resolution segmentation.<br>
<br>

## Sourcemodel creation
- EEG sourcemodel estimation by using CATs [[5]](#ref5) accurate cortex segmentation
<img align="right" width="250" src="images/sourcemodel_cat4000.png"> <br>
 and cortical thickness determination.
- A sourcemodel surface mesh is created at 2/3 distance between the gray and white matter surface using the cortical thickness estimation.
- The mesh vertices constitute the sourcemodel positions. The equivalent current dipoles are orientated according to their outwards pointing mesh normals.
<br>

## HArtMuT's artefact sourcemodel
- EEG ocular and muscular sourcemodel creation à la HArtMuT [[7]](#ref7)
<img align="right" width="250" src="images/hartmut.png"> <br>
is enabled by the additional creation of neck-extended meshes. 
- The small HArtMuT artefact sourcemodel is warped to the individual skull and scalp anatomy using a MATLAB translation of HArtMuT's [Individual Warp](https://github.com/harmening/HArtMuT/tree/main/individualwarp).
- The derived artefact sourcemodel allows not only computing cortical leadfields, but also ocular and muscular leadfields (that constitute major confounders in EEG data).
<br>

## Transformation to fiducial-based coordinate systems
 - Transformation of meshes, electrode positions and sourcemodels to coordinate systems, that are based on individual fiducical positions (e.g. CTF).
 - The individual fiducials (NAS, LPA, RPA) are gained through SPMs nonlinear mapping from the eTPM template.<br>
<br>

## Representing the cortical surface in Spherical Harmonics
 - The weighed spherical harmonic (weighted-SPHARM) representation can express the cortical surface as a weighted linear combination of spherical harmonics.
 - Since SPHARM is a better basis for the proper description of the cortex folding as the euclidean space, it is used in many different applications [[8]](#ref8).<br>
<br>


## Parcellation labels
- Cortical parcellation of the CAT12 cortical surface based registration
<img align="right" width="250" src="images/parcellation.png"> <br>
  of the hemsipheres to spheres. Examplarily integretated is the <br>
  Schaefer2018 atlas [[10]](#ref10).
- It is a widely used parcellation of the human cerebral cortex based on <br>
  intrinsic functional connectivity MRI.<br>
<br>



## Troubleshooting
- SPM12 sometimes fails due to fieldtrip routines replacing standard matlab
  routines. Make sure to delete (or remove from matlab search path) older matlab version folders from `./fieldtrip/compat/` and `./fieldtrip/external/spm12/external/fieldtrip/compat/`. For more information visit [https://github.com/fieldtrip/fieldtrip/blob/master/compat](https://github.com/fieldtrip/fieldtrip/blob/master/compat/matlablt2010b/README).
- Atra sometimes fails to correctly detect anatomical landmarks and one
  needs to manually specify the voxel coordinates of AC, PC and the vertex of
  the superior pontine sulcus (VSPS). Delete all files from your MRI folder
  exept the original MRI and the file `imagelist`. Create a new file called
  `landmarks.lm` and write the voxel coordinates of AC, PC and VSPS in the first, second and third line, respectively, seperated by one whitespace. Open `imagelist`, write the fullpath into the second line of `landmarks.lm` and restart the segmentation. For more details refer to [getting started with atra](https://www.nitrc.org/account/login.php?return_to=%2Ffrs%2Fdownload.php%2F10497%2Fgettingstartedwithatra1.0.pdf&feedback=The+tool%2Fresource+administrator+has+requested+that+you+log+in+to+download+this+file.) or [atra youtube tutorial](https://www.youtube.com/watch?v=q5GBaNnjOa8).
- `mri_convert` can fail on float-typed inputs (some MRIs, e.g. a few OASIS subjects, are distributed as float 4dfp) with `can only handle image types short and unsigned char`. Recast the volume to int16 before conversion/reorientation, e.g. with nibabel:
  ```python
  import nibabel as nib, numpy as np
  img = nib.load(fname)
  hdr = img.header
  hdr.set_data_dtype(np.int16)
  data = np.nan_to_num(img.get_fdata()).astype(np.int16)
  nib.save(nib.Nifti1Image(data, img.affine, hdr), fname)
  ```
<!--- 
mri = ft_read_mri('/fullpath/to/example.nii');
cfg = []; cfg.anaparameter = 'anatomy';
ft_sourceplot(cfg, mri);
--->
<br>
<br>


## Citing
If you find MRIsegmentation useful for your research, please consider citing:
```
@software{Harmening_MRIsegmentation_2022,
        author = {Harmening, Nils and Miklody, Daniel},
        doi = {10.5281/zenodo.7357674},
        month = {11},
        title = {{MRIsegmentation}},
        url = {https://github.com/harmening/MRIsegmentation},
        version = {1.1},
        year = {2022}
}
```
<br>
<br>

<!---
## Run full pipeline
The whole robust segmentation pipeline can be run via `main.m` or parallelized via the python wrapper `parallel.py`.<br>
<br>
<br>


## To Dos
- Make unit choosable
- Rewrite fieldtrip functions to kick out fieldtrip dependency
<br>
<br>
--->


<!---
## All matlab scripts + preprocessing can be run via main.m or parallelized via parallel.py
### It consists of
 - full segmentation pipeline call
 - surface boundary file conversion from matlab to tri (openmeeg/python)
(See also parallel_segmentation.py from Repo https://gitlab.tubit.tu-berlin.de/promillenille/OASIS)
**All matlab scripts can be run via the Huang_et_al_2013/run_segmentation_pipe.sh script. It consecutively calls**
### Huang_et_al_2013/start_seg.m
 - Gateway script to start running the SPM8 new segment.
 - Needed: SPM8 batch editor (https://www.fil.ion.ucl.ac.uk/spm/software/spm8/)
 - Input: T1, (T2), eTPM.nii
 - Output: Segmented tissue classes c1, c2, ... c6 into basefolder

### Huang_et_al_2013/mysegment.m
 - Automated Segmentation Clean Up + Tissue Smoothing (Andy Tools)
 - Needed: NIfTI Toolbox (http://www.mathworks.com/matlabcentral/fileexchange/8797)
 - Input: basefolder with c1-c6 files from SPM8
 - Output: Tissue masks mask_gray.nii, ... mask_skin.nii into basefolder

### MRIsegmentation/create_surface_meshes.m
 - Automatic triangle surface mesh generation. Tissues: scalp, skull, csf, cortex. Resolution each: 1922
 - Needed: FieldTrip Toolbox (http://www.fieldtriptoolbox.org/)
 - Input: basefolder with tissue masks from Andy tools
 - Output final segmentedmri.mat, surface meshes bnd4_1922.mat (in mm), bnd4_1922_corrected.mat (in m)

### MRIsegmentation/create_volume_meshes.m
 - Automatic volumetric tetrahedral mesh generation. Tissues: scalp, skull, csf, cortex. Voxel resolution: 1mm.
 - Needed: FieldTrip Toolbox (http://www.fieldtriptoolbox.org/)
 - Input: basefolder with tissue masks from Andy tools
 - Output final segmentedmri.mat, surface meshes bnd4_1922.mat (in mm), bnd4_1922_corrected.mat (in m)

### Huang_et_al_2013/electrode_placement.m (optional)
 - Electrode Placement of electrode coordinates D_74.mat and D_93.mat obtained from EasyCap (http://www.easycap.de/easycap/e/downloads/M1_ThetaPhi.htm)
 - Needed: NIfTI Toolbox (http://www.mathworks.com/matlabcentral/fileexchange/8797)<br />
           Natural Cubic Spline (http://www.mathworks.co.jp/matlabcentral/fileexchange/7617)
 - Input: D_74, D_93 or individual electrode coordinates, scalp_mask, nasion,inion,right,left,front_neck (required landmarks provided by the user through MRIcro 3D viewer)
 - Output: mask_elec.nii and mask_gel.nii 
 - CAUTION: not tested !
--->



<!---
### Not used: 
 - nipype
 - afni (https://afni.nimh.nih.gov/pub/dist/doc/htmldoc/background_install/install_instructs/steps_mac.html)
 - some more small scripts
--->


## References
<a id="ref1">[1]</a> Babak A. Ardekani and Alvin H. Bachman. "Model-based automatic detection of the anterior and posterior commissures on MRI scans" *NeuroImage*, (2009): **46**(3), 677 - 682. [``doi:10.1016/j.neuroimage.2009.02.030``](https://doi.org/10.1016/j.neuroimage.2009.02.030). <br> 
<!---  \cite{Ardekani09}  --->
<a id="ref2">[2]</a>  John Ashburner, Karl J. Friston. "Unified segmentation" *NeuroImage*, (2005): **26**(3), 839 - 851. [``doi:10.1016/j.neuroimage.2005.02.018``](https://doi.org/10.1016/j.neuroimage.2005.02.018). <br> 
<!--- \cite{Ashburner05} --->
<a id="ref3">[3]</a>  Yu Huang, Jacek Dmochowski, Yuzhuo Su, Abhishek Datta, Chris Rorden, Lucas Parra. "Automated MRI Segmentation for Individualized Modeling of Current Flow in the Human Head" *Journal of Neural Engineering*, (2013): **10**, 066004. [``doi:10.1088/1741-2560/10/6/066004``](https://doi.org/10.1088/1741-2560/10/6/066004). <br> 
<!--- \cite{Huang2013}  --->
<a id="ref4">[4]</a>  Qianqian Fang, David A. Boas. "Tetrahedral Mesh Generation from Volumetric Binary and Gray-Scale Images" *Proceedings of IEEE International Symposium on Biomedical Imaging*, (2009), 1142-1145. [``doi:10.1109/ISBI.2009.5193259``](https://doi.org/10.1109/ISBI.2009.5193259).<br> 
<!---  \cite{Fang2009*} --->
<a id="ref5">[5]</a> Christian Gaser, Robert Dahnke. "CAT - A Computational Anatomy Toolbox for the Analysis of Structural MRI Data" (2016): [http://www.neuro.uni-jena.de/hbm2016/GaserHBM2016.pdf](http://www.neuro.uni-jena.de/hbm2016/GaserHBM2016.pdf).<br>
<!---  \cite{Gaser2016} --->
<a id="ref6">[6]</a>  Robert Oostenveld, Pascal Fries, Eric Maris, Jan-Mathijs Schoffelen. "FieldTrip: Open Source Software for Advanced Analysis of MEG, EEG, and Invasive Electrophysiological Data" *Computational intelligence and neuroscience*, (2011), 156869. [``doi:10.1155/2011/156869``](https://doi.org/10.1155/2011/156869). <br> 
<!---  \cite{Fieldtrip11}  --->
<a id="ref7">[7]</a>  Nils Harmening, Marius Klug, Klaus Gramann, Daniel Miklody. "HArtMuT - modeling eye and muscle contributors in neuroelectric imaging" *Journal of Neural Engineering*, (2022), **19**066041. [``doi:10.1088/1741-2552/aca8ce``](https://doi.org/10.1088/1741-2552/aca8ce). <br> 
<!---  \cite{Harmening22}  --->
<a id="ref8">[8]</a>  Moo K. Chung, Kim M. Dalton, Li Shen, Alan C. Evans, Richard J. Davidson. "Weighted Fourier Series Representation and Its Application to Quantifying the Amount of Gray Matter" *IEEE Transactions on Mecical Imaging*, (2007): **26**(4). [``10.1109/TMI.2007.892519``](https://doi.org/10.1109/TMI.2007.892519). <br>
<!---  \cite{Chung2007}  --->
<a id="ref9">[9]</a>  Daniel S. Marcus, Tracy H. Wang, Jamie Parker, John G. Csernansky, John C. Morris, Randy L. Buckner. "Open Access Series of Imaging Studies (OASIS): Cross-sectional MRI Data in Young, Middle Aged, Nondemented, and Demented Older Adults" *Journal of Cognitive Neuroscience*, (2007): **19**(9), 1498-1507. [``doi:10.1162/jocn.2007.19.9.1498``](https://doi.org/10.1162/jocn.2007.19.9.1498). <br> 
<!--- \cite{OASIS2007}  --->
<a id="ref10">[10]</a>  Alexander Schaefer, Ru Kong, Evan M. Gordon, Timothy O. Laumann, Xi-Nian Zuo, Avram J. Holmes, Simon B. Eickhoff, B. T. Thomas Yeo. "Local-Global Parcellation of the Human Cerebral Cortex from Intrinsic Functional Connectivity MRI" *Cerebral Cortex*, (2018): **28**(9), 3095–3114. [``doi:10.1093/cercor/bhx179``](https://doi.org/10.1093/cercor/bhx179). <br>
<!---  \cite{Schaefer2018}  --->


