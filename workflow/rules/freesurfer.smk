"""FreeSurfer path: recon-all + mkheadsurf, exported into the NIfTI voxel grid.

Included by workflow/Snakefile when config freesurfer.enabled; uses its globals.
Runs in its own container image (config freesurfer.image), in parallel to the
MATLAB steps. Layout:
    results/<subject>/fs/<BASE>/          FreeSurfer subject dir (SUBJECTS_DIR = fs/)
    results/<subject>/fs/*.obj            surfaces in 1-based voxel indices of the
                                          staged NIfTI (see scripts/fs_surfaces_to_voxel.py)
"""

import struct
import sys

FS = config["freesurfer"]
FS_SD = f"{D}/fs"  # SUBJECTS_DIR
FS_SUBJ = f"{FS_SD}/{BASE}"
FS_EXPORTS = ["scalp", "lh.white", "rh.white", "lh.pial", "rh.pial", "cortex_white", "cortex_pial"]
FS_SMOOTH = int(FS.get("scalp_smooth_iterations", 20))
if FS_SMOOTH > 0:
    FS_EXPORTS.append("scalp_smoothed")

DEC = FS.get("decimate") or {}
DEC_ON = bool(DEC.get("enabled", False))
DEC_TARGETS = DEC.get("targets") or {}
DEC_SUPPORT = set(DEC.get("ribbon_support") or [])
# decimated surfaces, as <source>_<target vertex count>
DEC_OUTS = [f"{FS_SD}/{src}_{int(sched[-1])}.obj" for src, sched in DEC_TARGETS.items()]
# fsaverage6 Schaefer assets shipped with the repo (see scripts/fs_parcels.py)
FS_ATLAS_DIR = os.path.join(os.path.dirname(workflow.basedir), "Schaefer2018_Parcellations")
FS_ATLAS = FS.get("parcel_atlas", "Schaefer17_600")
# Decimated surfaces that carry a parcel per vertex: they need the voxel support
# feedback and annotations have to be configured. Only these can feed the head
# model export; the others (the scalp) are decimated plainly.
DEC_PARCELS = [s for s in DEC_TARGETS if s in DEC_SUPPORT] if DEC.get("annot") else []
DEC_PLAIN = [s for s in DEC_TARGETS if s not in DEC_PARCELS]

# cedalion head model export: needs the decimated cortex surface with its parcels
# and a decimated scalp surface (the smoothed one when it is being built)
HM = FS.get("headmodel") or {}
HM_SCALP = next((s for s in ("scalp_smoothed", "scalp") if s in DEC_TARGETS), None)
HM_ON = bool(HM.get("enabled", False)) and "cortex_pial" in DEC_PARCELS and HM_SCALP is not None
HM_DIR = f"{D}/headmodel"
HM_MASKS = ["mask_csf.nii", "mask_gray.nii", "mask_skin.nii", "mask_bone.nii", "mask_white.nii"]
# what cedalion.data.get_*_headmodel_files() reads, minus the two files that are
# still missing: cortex_pial_inflated.obj and landmarks.mrk.json
HM_OUTS = (
    [f"{HM_DIR}/{m}" for m in HM_MASKS]
    + [
        f"{HM_DIR}/cortex_pial_high.obj",
        f"{HM_DIR}/mask_brain.obj",
        f"{HM_DIR}/cortex_pial_inflated.obj",
        f"{HM_DIR}/mask_scalp.obj",
        f"{HM_DIR}/brain_vertex_coordinates.csv",
        f"{HM_DIR}/voxel_to_vertex_brain.mtx.gz",
        f"{HM_DIR}/parcel_colors.json",
    ]
    if HM_ON
    else []
)

# landmarks are built on the head model's own scalp, so they come after it; the
# fiducials they are constructed from are provided per subject (sample sheet)
LM = FS.get("landmarks") or {}
LM_ON = HM_ON and bool(LM.get("enabled", False))
LM_OUTS = [f"{HM_DIR}/landmarks.mrk.json"] if LM_ON else []


def lm_fiducials(wc):
    """The subject's fiducial markup file, from the sample sheet."""
    path = SAMPLES[wc.subject].get("landmarks")
    if not path:
        raise ValueError(
            f"{wc.subject}: freesurfer.landmarks is enabled but the sample sheet has no "
            "'landmarks' column entry (a .mrk.json with Nz, Iz, LPA, RPA; see "
            "config/landmarks/README.md)"
        )
    return path


def hm_dec(surface):
    """Path of a decimated surface at its final vertex count."""
    return f"{FS_SD}/{surface}_{int(DEC_TARGETS[surface][-1])}.obj"


def hm_mni_from(wc):
    """config headmodel.mni_from, overridable per subject in the sample sheet."""
    return SAMPLES[wc.subject].get("mni_from") or HM.get("mni_from", "talairach")
DEC_SCRIPT = os.path.join(workflow.basedir, "scripts", "fs_decimate_mesh.py")
DEC_PYTHON = DEC.get("python") or sys.executable


def dec_schedule(surface):
    return DEC_TARGETS[surface]


def dec_extra_args(wc, input):
    """Voxel-support and annotation arguments shared by both decimation rules."""
    args = ""
    if input.ribbon:
        args += (
            f" --reference {input.reference[0]} --ribbon {input.ribbon[0]}"
            f" --voxel-stealing-from {int(DEC.get('voxel_stealing_from', 26000))}"
            f" --quality-radius {DEC.get('quality_radius', 15)}"
        )
    if input.annots:
        args += f" --annot-lh {input.annots[0]} --annot-rh {input.annots[1]}"
    return args


def dec_annots(wc):
    """Surface annotation files for parcel-aware decimation (empty if unused)."""
    pattern = DEC.get("annot")
    if not pattern or wc.surface not in DEC_SUPPORT:
        return []
    subj = FS_SUBJ.format(subject=wc.subject)
    return [pattern.format(fs_subj=subj, hemi=h) for h in ("lh", "rh")]


def dec_ribbon(wc):
    """ribbon.mgz, only for the surfaces that need voxel support."""
    if wc.surface not in DEC_SUPPORT:
        return []
    return [f"{FS_SUBJ}/mri/ribbon.mgz".format(subject=wc.subject)]


def dec_reference(wc):
    return [IMG.format(subject=wc.subject)] if wc.surface in DEC_SUPPORT else []


def fs_runner():
    lic = FS.get("license") or ""
    if not lic:
        raise ValueError("config freesurfer.license: path to the FreeSurfer license.txt is required")
    scripts = os.path.join(workflow.basedir, "scripts")
    binds = [(A(RES), A(RES)), (A(lic), A(lic)), (scripts, scripts)]
    env = [("FS_LICENSE", A(lic))]
    return runner(FS["image"], binds, env)


# recon-all & co. rely on the environment SetUpFreeSurfer.sh exports
FSRUN = fs_runner() + ' bash -c \'source "$FREESURFER_HOME/SetUpFreeSurfer.sh" >/dev/null 2>&1; '
FS_SCRIPT = os.path.join(workflow.basedir, "scripts", "fs_surfaces_to_voxel.py")


def nifti_fov_mm(path):
    """Largest field of view (mm) of a NIfTI-1 image, from the header only."""
    with open(path, "rb") as fh:
        hdr = fh.read(348)
    endian = "<" if struct.unpack("<i", hdr[:4])[0] == 348 else ">"
    dim = struct.unpack(endian + "8h", hdr[40:56])
    pixdim = struct.unpack(endian + "8f", hdr[76:108])
    return max(abs(dim[i] * pixdim[i]) for i in range(1, 4))


rule recon_all:
    """FreeSurfer recon-all (all stages), T2 for pial refinement when available."""
    input:
        img=IMG,
        t2=t2_input(IMG_T2),
    output:
        done=f"{FS_SUBJ}/scripts/recon-all.done",
        t1=f"{FS_SUBJ}/mri/T1.mgz",
        surfs=expand(f"{FS_SUBJ}/surf/{{h}}.{{s}}", h=["lh", "rh"], s=["white", "pial", "sphere.reg"], allow_missing=True),
    log:
        f"{LOGS}/recon_all.log",
    params:
        sd=lambda wc: A(FS_SD.format(subject=wc.subject)),
        subj=lambda wc: BASE.format(subject=wc.subject),
        t2=lambda wc, input: f"-T2 {A(input.t2[0])} -T2pial" if input.t2 and FS.get("use_t2", True) else "",
        # recon-all refuses fields of view above 256 mm unless told to crop
        cw256=lambda wc, input: "-cw256" if nifti_fov_mm(input.img) > 256 else "",
        extra=FS.get("recon_all_args", "") or "",
    threads: int(FS.get("threads", 8))
    resources:
        mem_mb=int(FS.get("mem_mb", 12000)),
    shell:
        # recon-all -i refuses to start if the subject directory already exists
        "rm -rf {params.sd}/{params.subj} && mkdir -p {params.sd} && "
        + FSRUN
        + "recon-all -sd {params.sd} -s {params.subj} -i $(realpath {input.img}) {params.t2} {params.cw256} "
        "-all -parallel -openmp {threads} {params.extra}' > {log} 2>&1"


rule head_surface:
    """Scalp surface from the conformed T1 (mkheadsurf), cleaned of thin artefact sheets, optionally smoothed."""
    input:
        done=f"{FS_SUBJ}/scripts/recon-all.done",
        t1=f"{FS_SUBJ}/mri/T1.mgz",
    output:
        f"{FS_SUBJ}/surf/lh.seghead",
        f"{FS_SUBJ}/surf/lh.smseghead" if FS_SMOOTH > 0 else [],
        mask=f"{FS_SUBJ}/mri/seghead.mgz",
        raw=f"{FS_SUBJ}/mri/seghead.orig.mgz",
    log:
        f"{LOGS}/head_surface.log",
    params:
        sd=lambda wc: A(FS_SD.format(subject=wc.subject)),
        subj=lambda wc: BASE.format(subject=wc.subject),
        clean=os.path.join(workflow.basedir, "scripts", "fs_clean_headmask.py"),
        radius=int(FS.get("scalp_open_radius", 1)),
        prior=lambda wc, output: (
            ""
            if FS.get("scalp_margin_mm") is None
            else " --talairach {0}/mri/transforms/talairach.xfm --t1 {0}/mri/T1.mgz --margin {1} --dim-below {2}".format(
                os.path.dirname(os.path.dirname(A(output.mask))),
                FS.get("scalp_margin_mm", 2),
                FS.get("scalp_dim_below", 70),
            )
        ),
        smooth=lambda wc, output: (
            f" && mris_smooth -nw -n {FS_SMOOTH} {A(output[0])} {A(output[1])}" if FS_SMOOTH > 0 else ""
        ),
    threads: 1
    resources:
        mem_mb=4000,
    shell:
        # mkheadsurf = mri_seghead + mri_mc + mris_smooth; the head mask is cleaned
        # in between and the surface re-tessellated the same way mkheadsurf does
        FSRUN
        + "export SUBJECTS_DIR={params.sd}; mkheadsurf -subjid {params.subj} -srcvol T1.mgz"
        " && cp {output.mask} {output.raw}"
        " && fspython {params.clean} {output.raw} {output.mask} --radius {params.radius}{params.prior}"
        " && mri_mc {output.mask} 1 {output[0]} && mris_smooth -n 10 -nw {output[0]} {output[0]}"
        "{params.smooth}' > {log} 2>&1"


rule fs_export:
    """Scalp and cortical surfaces as OBJ in 1-based voxel indices of the staged NIfTI."""
    input:
        img=IMG,
        t1=f"{FS_SUBJ}/mri/T1.mgz",
        surfs=expand(f"{FS_SUBJ}/surf/{{h}}.{{s}}", h=["lh", "rh"], s=["white", "pial"], allow_missing=True),
        scalp=f"{FS_SUBJ}/surf/lh.seghead",
        scalp_smoothed=f"{FS_SUBJ}/surf/lh.smseghead" if FS_SMOOTH > 0 else [],
    output:
        expand(f"{FS_SD}/{{name}}.obj", name=FS_EXPORTS, allow_missing=True),
    log:
        f"{LOGS}/fs_export.log",
    params:
        script=FS_SCRIPT,
        out=lambda wc: A(FS_SD.format(subject=wc.subject)),
        specs=lambda wc, input: " ".join(
            [f"scalp={A(input.scalp)}"]
            + ([f"scalp_smoothed={A(str(input.scalp_smoothed))}"] if input.scalp_smoothed else [])
            + [f"{os.path.basename(s)}={A(s)}" for s in input.surfs]
        ),
    threads: 1
    resources:
        mem_mb=2000,
    shell:
        # fspython ships numpy + nibabel, so the export runs in the same image
        fs_runner()
        + " fspython {params.script} $(realpath {input.img}) $(realpath {input.t1}) {params.out} {params.specs} > {log} 2>&1"


rule compare_surfaces:
    """QC: exact cross-sections + distances of MATLAB vs FreeSurfer scalp and cortex."""
    input:
        img=IMG,
        m_scalp=f"{D}/mask_scalp.obj",
        m_cortex=f"{D}/mask_brain.obj",
        f_scalp=f"{FS_SD}/scalp{'_smoothed' if FS_SMOOTH > 0 else ''}.obj",
        f_pial=f"{FS_SD}/cortex_pial.obj",
        f_white=f"{FS_SD}/cortex_white.obj",
    output:
        png=f"{RES}/qc/{{subject}}/surfaces_matlab_vs_freesurfer.png",
        json=f"{RES}/qc/{{subject}}/surfaces_matlab_vs_freesurfer.json",
    log:
        f"{LOGS}/compare_surfaces.log",
    params:
        script=os.path.join(workflow.basedir, "scripts", "compare_surfaces.py"),
        python=sys.executable,
    localrule: True
    shell:
        # runs with Snakemake's own interpreter: needs numpy, scipy, nibabel, matplotlib
        "{params.python} {params.script} {input.img} {output.png} {output.json} "
        "'MATLAB scalp=scalp:lime:{input.m_scalp}' 'MATLAB cortex (CAT12 central)=cortex:magenta:{input.m_cortex}' "
        "'FreeSurfer scalp=scalp:cyan:{input.f_scalp}' 'FreeSurfer pial=cortex:red:{input.f_pial}' "
        "'FreeSurfer white=cortex:yellow:{input.f_white}' > {log} 2>&1"


rule compare_masks:
    """QC: volumes, Dice and cross-tabulation of the MATLAB tissue masks vs FreeSurfer's aseg."""
    input:
        reference=IMG,
        masks=MASKS,
        aseg=f"{FS_SUBJ}/mri/aseg.mgz",
        head=f"{FS_SUBJ}/mri/seghead.mgz",
    output:
        png=f"{RES}/qc/{{subject}}/masks_matlab_vs_freesurfer.png",
        json=f"{RES}/qc/{{subject}}/masks_matlab_vs_freesurfer.json",
    log:
        f"{LOGS}/compare_masks.log",
    params:
        script=os.path.join(workflow.basedir, "scripts", "compare_masks.py"),
        python=sys.executable,
        masks=lambda wc: D.format(subject=wc.subject),
    localrule: True
    shell:
        # runs with Snakemake's own interpreter: needs numpy, nibabel, matplotlib
        "{params.python} {params.script} --reference {input.reference} --masks {params.masks}"
        " --aseg {input.aseg} --head {input.head}"
        " --out-png {output.png} --out-json {output.json} > {log} 2>&1"


if DEC_ON and "cortex_pial" in DEC_TARGETS:

    rule compare_voxel_grouping:
        """QC: SPM grey matter vs FreeSurfer ribbon as the voxel set for the cortex surface."""
        input:
            reference=IMG,
            surface=f"{FS_SD}/cortex_pial_{int(DEC_TARGETS['cortex_pial'][-1])}.obj",
            high_res=f"{FS_SD}/cortex_pial.obj",
            grey=f"{D}/mask_gray.nii",
            aseg=f"{FS_SUBJ}/mri/aseg.mgz",
            ribbon=f"{FS_SUBJ}/mri/ribbon.mgz",
            spheres=expand(f"{FS_SUBJ}/surf/{{h}}.sphere.reg", h=["lh", "rh"], allow_missing=True),
            atlas=FS_ATLAS_DIR,
        output:
            png=f"{RES}/qc/{{subject}}/voxel_grouping_matlab_vs_freesurfer.png",
            json=f"{RES}/qc/{{subject}}/voxel_grouping_matlab_vs_freesurfer.json",
        log:
            f"{LOGS}/compare_voxel_grouping.log",
        params:
            script=os.path.join(workflow.basedir, "scripts", "compare_voxel_grouping.py"),
            python=DEC.get("python") or sys.executable,
            fs_subject=lambda wc: FS_SUBJ.format(subject=wc.subject),
            atlas=FS_ATLAS,
        threads: 1
        resources:
            mem_mb=int(DEC.get("mem_mb", 16000)),
        shell:
            # needs cedalion, like fs_decimate
            "{params.python} {params.script} --reference {input.reference} --surface {input.surface}"
            " --fs-subject {params.fs_subject} --grey {input.grey} --atlas-dir {input.atlas}"
            " --atlas {params.atlas} --out-png {output.png} --out-json {output.json} > {log} 2>&1"


if DEC_ON and DEC.get("annot"):

    rule fs_annot:
        """Transfer the fsaverage6 Schaefer parcellation onto the subject's surface.

        Produces the annotations the parcel-aware decimation and the head model
        export read; recon-all itself only ships aparc/BA parcellations.
        """
        input:
            spheres=expand(f"{FS_SUBJ}/surf/{{h}}.sphere.reg", h=["lh", "rh"], allow_missing=True),
            atlas=FS_ATLAS_DIR,
        output:
            annots=[DEC["annot"].format(fs_subj=FS_SUBJ, hemi=h) for h in ("lh", "rh")],
        log:
            f"{LOGS}/fs_annot.log",
        params:
            script=os.path.join(workflow.basedir, "scripts", "fs_parcels.py"),
            python=DEC.get("python") or sys.executable,
            fs_subject=lambda wc: FS_SUBJ.format(subject=wc.subject),
            atlas=FS_ATLAS,
        threads: 1
        resources:
            mem_mb=int(DEC.get("mem_mb", 16000)),
        shell:
            # runs on the host like fs_decimate, but needs only nibabel and scipy
            "{params.python} {params.script} --fs-subject {params.fs_subject}"
            " --atlas-dir {input.atlas} --atlas {params.atlas}"
            " --out-lh {output.annots[0]} --out-rh {output.annots[1]} > {log} 2>&1"


if DEC_PLAIN:

    rule fs_decimate:
        """Decimate a FreeSurfer surface without the voxel feedback (e.g. the scalp)."""
        input:
            mesh=f"{FS_SD}/{{surface}}.obj",
            reference=dec_reference,
            ribbon=dec_ribbon,
            annots=dec_annots,
        output:
            mesh=f"{FS_SD}/{{surface}}_{{n,[0-9]+}}.obj",
        log:
            f"{LOGS}/fs_decimate_{{surface}}_{{n}}.log",
        wildcard_constraints:
            surface="|".join(re.escape(s) for s in DEC_PLAIN) or "(?!)",
        params:
            script=DEC_SCRIPT,
            python=DEC_PYTHON,
            nvertices=lambda wc: " ".join(str(int(n)) for n in dec_schedule(wc.surface)),
            extra=dec_extra_args,
        threads: 1
        resources:
            mem_mb=int(DEC.get("mem_mb", 16000)),
        shell:
            "{params.python} {params.script} --mesh {input.mesh} --out-mesh {output.mesh}"
            " --nvertices {params.nvertices}{params.extra} > {log} 2>&1"


if DEC_PARCELS:

    rule fs_decimate_parcels:
        """Decimate a cortex surface, keeping each vertex's voxel support and parcel.

        Also keeps the voxel-to-vertex mapping of the last step, which the head
        model export turns into cedalion's voxel_to_vertex_brain matrix.
        """
        input:
            mesh=f"{FS_SD}/{{surface}}.obj",
            reference=dec_reference,
            ribbon=dec_ribbon,
            annots=dec_annots,
        output:
            mesh=f"{FS_SD}/{{surface}}_{{n,[0-9]+}}.obj",
            parcels=f"{FS_SD}/{{surface}}_{{n,[0-9]+}}_parcels.json",
            vertexinfo=f"{FS_SD}/{{surface}}_{{n,[0-9]+}}_vertexinfo.npz",
        log:
            f"{LOGS}/fs_decimate_{{surface}}_{{n}}.log",
        wildcard_constraints:
            surface="|".join(re.escape(s) for s in DEC_PARCELS) or "(?!)",
        params:
            script=DEC_SCRIPT,
            python=DEC_PYTHON,
            nvertices=lambda wc: " ".join(str(int(n)) for n in dec_schedule(wc.surface)),
            extra=dec_extra_args,
        threads: 1
        resources:
            mem_mb=int(DEC.get("mem_mb", 16000)),
        shell:
            "{params.python} {params.script} --mesh {input.mesh} --out-mesh {output.mesh}"
            " --nvertices {params.nvertices}{params.extra}"
            " --out-parcels {output.parcels} --out-vertexinfo {output.vertexinfo} > {log} 2>&1"


if HM_ON:

    rule headmodel:
        """Assemble results/<subject>/headmodel/, a cedalion head model directory."""
        input:
            reference=IMG,
            high_res=f"{FS_SD}/cortex_pial.obj",
            reduced=hm_dec("cortex_pial"),
            vertexinfo=hm_dec("cortex_pial").replace(".obj", "_vertexinfo.npz"),
            scalp=hm_dec(HM_SCALP),
            inflated=expand(f"{FS_SUBJ}/surf/{{h}}.inflated", h=["lh", "rh"], allow_missing=True),
            masks=[f"{D}/{m}" for m in HM_MASKS],
            talairach=f"{FS_SUBJ}/mri/transforms/talairach.xfm",
            atlas=FS_ATLAS_DIR,
        output:
            HM_OUTS,
        log:
            f"{LOGS}/headmodel.log",
        params:
            script=os.path.join(workflow.basedir, "scripts", "make_headmodel.py"),
            python=DEC_PYTHON,
            out_dir=HM_DIR,
            atlas=FS_ATLAS,
            mni_from=hm_mni_from,
        threads: 1
        resources:
            mem_mb=int(DEC.get("mem_mb", 16000)),
        shell:
            # needs cedalion, like fs_decimate
            "{params.python} {params.script} --out-dir {params.out_dir}"
            " --reference {input.reference} --high-res {input.high_res}"
            " --reduced {input.reduced} --vertexinfo {input.vertexinfo}"
            " --scalp {input.scalp} --inflated {input.inflated} --masks {input.masks}"
            " --atlas-dir {input.atlas} --atlas {params.atlas}"
            " --mni-from {params.mni_from} --talairach {input.talairach} > {log} 2>&1"


if LM_ON:

    rule landmarks:
        """Build the 10-10 landmark system on the head model's scalp surface."""
        input:
            scalp=f"{HM_DIR}/mask_scalp.obj",
            reference=IMG,
            fiducials=lm_fiducials,
        output:
            f"{HM_DIR}/landmarks.mrk.json",
        log:
            f"{LOGS}/landmarks.log",
        params:
            script=os.path.join(workflow.basedir, "scripts", "make_landmarks.py"),
            python=DEC_PYTHON,
            snap=" ".join(LM.get("snap", ["Nz", "Iz", "LPA", "RPA"])),
            max_snap=float(LM.get("max_snap_mm", 20)),
        threads: 1
        resources:
            mem_mb=int(DEC.get("mem_mb", 16000)),
        shell:
            # needs cedalion, like fs_decimate
            "{params.python} {params.script} --scalp {input.scalp}"
            " --reference {input.reference} --fiducials {input.fiducials}"
            " --out {output} --max-snap-mm {params.max_snap}"
            " --snap {params.snap} > {log} 2>&1"

    rule compare_landmarks:
        """QC: the constructed landmarks on the scalp, against a reference set."""
        input:
            landmarks=f"{HM_DIR}/landmarks.mrk.json",
            scalp=f"{HM_DIR}/mask_scalp.obj",
            reference=IMG,
        output:
            png=f"{RES}/qc/{{subject}}/landmarks.png",
            json=f"{RES}/qc/{{subject}}/landmarks.json",
        log:
            f"{LOGS}/compare_landmarks.log",
        params:
            script=os.path.join(workflow.basedir, "scripts", "compare_landmarks.py"),
            python=DEC_PYTHON,
            # 'fieldtrip' downloads the reference through cedalion rather than
            # shipping it with the workflow; empty = no comparison
            ref=lambda wc: SAMPLES[wc.subject].get("landmarks_reference") or "",
        threads: 1
        resources:
            mem_mb=int(DEC.get("mem_mb", 16000)),
        shell:
            "{params.python} {params.script} --landmarks {input.landmarks}"
            " --scalp {input.scalp} --reference {input.reference} --subject {wildcards.subject}"
            " --reference-landmarks '{params.ref}'"
            " --out-png {output.png} --out-json {output.json} > {log} 2>&1"
