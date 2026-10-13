"""EEG pipeline rules (start_segmentation.m): electrodes, HArtMuT, sourcemodel, CTF, SPHARM.

Included by workflow/Snakefile when config pipeline == "eeg"; uses its globals.
"""


rule align_electrodes:
    """Nonlinearly warp template electrodes/fiducials and project them onto the scalp."""
    input:
        img=IMG,
        bnd=f"{D}/bnd4_{NV}_corrected.mat",
        warp=f"{D}/y_{BASE}.nii",
    output:
        expand(f"{D}/elec_standard_{{m}}.mat", m=["1005", "1010", "1020"], allow_missing=True),
        f"{D}/fiducials_nonlinear.mat",
    log:
        f"{LOGS}/align_electrodes.log",
    params:
        cmd=lambda wc, input: f"S = load('{A(input.bnd)}'); align_nonlinear('{A(input.img)}', S.new_bnd(1))",
    threads: THREADS
    resources:
        mem_mb=MEM_MB,
    shell:
        mriseg("{params.cmd}")


rule hartmut_sources:
    """Warp HArtMuT's ocular/muscular artefact sources to the individual head (ACPC only)."""
    input:
        scalp=f"{D}/bnd4_{NV}_hartmut_corrected_scalp.tri",
        skull=f"{D}/bnd4_{NV}_hartmut_corrected_skull.tri",
    output:
        f"{D}/artefact_sourcemodel_HArtMuT_small.mat",
    log:
        f"{LOGS}/hartmut_sources.log",
    params:
        cmd=lambda wc, input, output: (
            f"h = {HARTMUT}; "
            "warp_HArtMuT_sources("
            "fullfile(h, 'individualwarp', 'NYhead', 'scalp.stl'), "
            "fullfile(h, 'individualwarp', 'NYhead', 'skull.stl'), "
            f"'{A(input.scalp)}', '{A(input.skull)}', "
            "fullfile(h, 'HArtMuTmodels', 'HArtMuT_NYhead_small.mat'), "
            f"'{A(output[0])}', [0 0 0])"
        ),
    threads: THREADS
    resources:
        mem_mb=MEM_MB,
    shell:
        mriseg("{params.cmd}")


rule hartmut_bnd:
    """Bundle neck-extended surfaces (+ warped eye sources) for the FieldTrip HArtMuT example."""
    input:
        img=IMG,
        bnd=f"{D}/bnd4_{NV}_hartmut_corrected.mat",
        sources=f"{D}/artefact_sourcemodel_HArtMuT_small.mat",
    output:
        f"{D}/mri_bnd.mat",
    log:
        f"{LOGS}/hartmut_bnd.log",
    params:
        cmd=lambda wc, input: f"export_fieldtrip_hartmut_bnd('{A(input.img)}', {NV}, 'mri_bnd.mat')",
    threads: THREADS
    resources:
        mem_mb=MEM_MB,
    shell:
        mriseg("{params.cmd}")


rule sourcemodel:
    """Cortical sourcemodel at 2/3 GM/WM depth from the CAT12 surfaces."""
    input:
        img=IMG,
        surfaces=CAT_SURFACES,
    output:
        f"{D}/sourcemodel{NS}.mat",
        f"{D}/sourcemodel{NS}.tri",
        f"{D}/mask_brain{NC}.stl",
        expand(f"{CAT_SURF}/{{h}}_{{s}}.mat", h=["lh", "rh"], s=["surf", "sphere"], allow_missing=True),
    log:
        f"{LOGS}/sourcemodel.log",
    params:
        cmd=lambda wc, input: f"read_cat({CAT12}, '{A(input.img)}', {NS}, {NC})",
    threads: THREADS
    resources:
        mem_mb=MEM_MB,
    shell:
        mriseg("{params.cmd}")


rule transform_to_ctf:
    """Transform meshes, electrodes and sourcemodels into the fiducial-based CTF system."""
    input:
        img=IMG,
        fiducials=f"{D}/fiducials_nonlinear.mat",
        elec=expand(f"{D}/elec_standard_{{m}}.mat", m=["1005", "1010", "1020"], allow_missing=True),
        bnd=f"{D}/bnd4_{NV}_corrected.mat",
        bnd_hartmut=f"{D}/bnd4_{NV}_hartmut_corrected.mat",
        mesh=f"{D}/mesh6_maxvox{MV}_hartmut.mat",
        sourcemodel=f"{D}/sourcemodel{NS}.mat",
        hemis=expand(f"{CAT_SURF}/{{h}}_{{s}}.mat", h=["lh", "rh"], s=["surf", "sphere"], allow_missing=True),
    output:
        f"{D}/ctf/acpc2ctf.mat",
        f"{D}/ctf/sourcemodel{NS}.mat",
        expand(f"{D}/ctf/elec_standard_{{m}}.mat", m=["1005", "1010", "1020"], allow_missing=True),
        expand(f"{D}/ctf/{{h}}_{{s}}.mat", h=["lh", "rh"], s=["surf", "sphere"], allow_missing=True),
    log:
        f"{LOGS}/transform_to_ctf.log",
    params:
        cmd=lambda wc, input: f"transform_to_ctf('{A(input.img)}', {NV}, {MV}, {NS})",
    threads: THREADS
    resources:
        mem_mb=MEM_MB,
    shell:
        mriseg("{params.cmd}")


rule spharm:
    """Weighted spherical harmonics representation of both hemispheres (CTF frame)."""
    input:
        img=IMG,
        hemis=expand(f"{D}/ctf/{{h}}_{{s}}.mat", h=["lh", "rh"], s=["surf", "sphere"], allow_missing=True),
    output:
        expand(f"{D}/ctf/{{h}}_fourier_{{k}}.mat", h=["lh", "rh"], k=["coeff", "vector"], allow_missing=True),
    log:
        f"{LOGS}/spharm.log",
    params:
        cmd=lambda wc, input: f"spharm('{A(input.img)}', {SPHARM_L})",
    threads: THREADS
    resources:
        mem_mb=MEM_MB,
    shell:
        mriseg("{params.cmd}")
