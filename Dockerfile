# MRIsegmentation container image
#
# Bundles MATLAB (installed with mpm), SPM12, CAT12, FieldTrip, iso2mesh,
# Andy's tools (Huang et al. 2013), the NIfTI toolbox, HArtMuT, weighted-SPHARM
# and this repository's pipeline code. Every pipeline step can be run
# non-interactively through the `mriseg` wrapper, which is what the Snakemake
# workflow in workflow/ does.
#
# Build (from the repository root; spm12.zip, cat12_r2159.zip and
# Huang_et_al_2013.zip must be present, see docker/fetch-deps.sh):
#   podman build -t mrisegmentation .
#   docker build -t mrisegmentation .
#
# Run one pipeline step with apptainer (online licensing via your MathWorks
# account, see below; `bash docker/build-sif.sh` converts the image to mriseg.sif):
#   apptainer exec mriseg.sif \
#       mriseg "start_seg('/path/to/subject/T1.nii', [], 'Huang_et_al_2013/eTPM.nii', false)"
# ... or with podman/docker and a network license / batch token. With rootless
# podman, `--user root` makes output files belong to you (container root == host
# user); with docker use `--user $(id -u):$(id -g)` instead:
#   podman run --rm --user root -e MLM_LICENSE_FILE \
#       -v /path/to/subject:/path/to/subject mrisegmentation mriseg "..."
#
# Interactive MATLAB shell:
#   apptainer exec mriseg.sif matlab -nodesktop
#
# Licensing: MATLAB inside the container needs one of
#   * your MathWorks account (Individual/online licensing, the default build):
#     run the image with apptainer as yourself, it reuses the sign-in cached by
#     the MATLAB on your host (see docker/build-sif.sh and the README)
#   * a batch licensing token         -> -e MLM_LICENSE_TOKEN=user@mail::token
#   * a network license manager       -> -e MLM_LICENSE_FILE=27000@server
#   * a license file                  -> -v lic.lic:/licenses/lic.lic -e MLM_LICENSE_FILE=/licenses/lic.lic
#   (the last three need --build-arg MATLAB_LICMODE=file)
# See https://github.com/mathworks-ref-arch/matlab-dockerfile/tree/main/alternates/non-interactive

# ============================================================================
# Build arguments
# ============================================================================
ARG MATLAB_RELEASE=R2025b
# Toolboxes needed by the pipeline. Image Processing is mandatory (Andy's tools,
# mesh creation); the others are used by FieldTrip/CAT12 code paths. Trim with
# --build-arg MATLAB_TOOLBOXES="Image_Processing_Toolbox" to get a smaller image.
ARG MATLAB_TOOLBOXES="Image_Processing_Toolbox Signal_Processing_Toolbox Statistics_and_Machine_Learning_Toolbox Optimization_Toolbox Parallel_Computing_Toolbox Wavelet_Toolbox"
ARG MATLAB_INSTALL_LOCATION=/opt/matlab/${MATLAB_RELEASE}

# ============================================================================
# Base image: MathWorks' dependency image (no MATLAB yet, installed via mpm)
# ============================================================================
FROM mathworks/matlab-deps:${MATLAB_RELEASE}

ARG MATLAB_RELEASE
ARG MATLAB_TOOLBOXES
ARG MATLAB_INSTALL_LOCATION

LABEL org.opencontainers.image.title="MRIsegmentation" \
      org.opencontainers.image.description="MRI segmentation pipeline (SPM12, CAT12, FieldTrip, iso2mesh) for FEM/BEM head models" \
      org.opencontainers.image.source="https://github.com/harmening/MRIsegmentation"

ENV DEBIAN_FRONTEND=noninteractive

# ============================================================================
# System dependencies
# ============================================================================
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
        ca-certificates \
        curl \
        git \
        patch \
        sudo \
        tini \
        unzip \
        wget \
        xvfb \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Unprivileged user that runs MATLAB
RUN adduser --shell /bin/bash --disabled-password --gecos "" matlab \
    && echo "matlab ALL=(ALL) NOPASSWD: ALL" > /etc/sudoers.d/matlab \
    && chmod 0440 /etc/sudoers.d/matlab

# ============================================================================
# MATLAB + toolboxes via MATLAB Package Manager (mpm). Largest layer, so it
# comes first and is only rebuilt when the release/toolbox list changes.
# ============================================================================
RUN wget -q https://www.mathworks.com/mpm/glnxa64/mpm -O /tmp/mpm \
    && chmod +x /tmp/mpm \
    && /tmp/mpm install \
        --release=${MATLAB_RELEASE} \
        --destination=${MATLAB_INSTALL_LOCATION} \
        --no-gpu \
        --products MATLAB ${MATLAB_TOOLBOXES} \
       || (echo "MPM installation failed:" && cat /tmp/mathworks_root.log && false) \
    && rm -f /tmp/mpm /tmp/mathworks_root.log \
    && ln -s ${MATLAB_INSTALL_LOCATION}/bin/matlab /usr/local/bin/matlab

# matlab-batch: runs MATLAB non-interactively with a batch licensing token
# (MLM_LICENSE_TOKEN). Installed to /opt/mathworks, symlinked to /usr/local/bin.
RUN wget -q https://raw.githubusercontent.com/mathworks-ref-arch/matlab-dockerfile/main/alternates/non-interactive/install/install-matlab-batch.sh \
        -O /tmp/install-matlab-batch.sh \
    && bash /tmp/install-matlab-batch.sh \
    && rm -f /tmp/install-matlab-batch.sh

# License mode. "onlinelicensing" writes the same marker a desktop install gets
# when you pick "Sign in with MathWorks account": MATLAB then licenses through
# your account's cached sign-in (~/.matlab, ~/.MathWorks of the calling user, so
# run the image with apptainer, which shares your home) and IGNORES
# MLM_LICENSE_FILE. Build with --build-arg MATLAB_LICMODE=file for a network
# license manager, license file or batch token instead.
ARG MATLAB_LICMODE=onlinelicensing
RUN if [ "${MATLAB_LICMODE}" = "onlinelicensing" ]; then \
        mkdir -p ${MATLAB_INSTALL_LOCATION}/licenses \
        && printf '<?xml version="1.0" encoding="UTF-8" standalone="no" ?>\n<root>\n  <ActivationEntry hostname="*" matlabroot="*" user="*">\n    <licmode>onlinelicensing</licmode>\n  </ActivationEntry>\n</root>\n' \
            > ${MATLAB_INSTALL_LOCATION}/licenses/license_info.xml; \
    fi

ENV MATLAB_INSTALL_LOCATION=${MATLAB_INSTALL_LOCATION}
ENV PATH=${MATLAB_INSTALL_LOCATION}/bin:${PATH}
# Shell for MATLAB's system() calls (CAT12 binaries, ATRA); see docker/mriseg
ENV MATLAB_SHELL=/bin/bash

# ============================================================================
# Third-party MATLAB toolboxes (in dependency order of the README)
# ============================================================================
ENV MRISEG_HOME=/opt/MRIsegmentation
WORKDIR ${MRISEG_HOME}

# FieldTrip (ships SPM12 in external/, which is replaced by a full SPM12 below)
RUN git clone --depth 1 https://github.com/fieldtrip/fieldtrip.git

# SPM12 (bind-mounted from the build context; download with docker/fetch-deps.sh)
RUN --mount=type=bind,source=spm12.zip,target=/tmp/spm12.zip \
    rm -rf ./fieldtrip/external/spm12 \
    && unzip -q /tmp/spm12.zip -d ./fieldtrip/external/

# Rebuild the SPM12 MEX files against the installed MATLAB release. SPM12 ships
# precompiled glnxa64 binaries, so a failed rebuild keeps those instead of
# breaking the image (no `make distclean`, which would delete them first).
RUN cd ./fieldtrip/external/spm12/src \
    && (make && make install && make external && make external-install) \
       || echo "WARNING: SPM12 MEX rebuild failed, keeping the precompiled binaries"

# CAT12 r2159 (version pinned for compatibility with the pipeline's CAT settings)
RUN --mount=type=bind,source=cat12_r2159.zip,target=/tmp/cat12_r2159.zip \
    rm -rf ./fieldtrip/external/spm12/toolbox/cat12 \
    && unzip -q /tmp/cat12_r2159.zip -d ./fieldtrip/external/spm12/toolbox/

# iso2mesh (full version; FieldTrip only carries a subset)
RUN git clone --depth 1 https://github.com/fangq/iso2mesh.git

# Andy's tools (Huang et al. 2013) incl. eTPM.nii, patched to store the
# nonlinear warp fields (needed for electrode/fiducial alignment) and to
# survive heads with few connected CSF/skin components (docker/patches/)
RUN --mount=type=bind,source=Huang_et_al_2013.zip,target=/tmp/Huang_et_al_2013.zip \
    --mount=type=bind,source=docker/patches,target=/tmp/patches \
    unzip -q /tmp/Huang_et_al_2013.zip -d ./Huang_et_al_2013 \
    && sed -i 's/warp\.write = \[0 0\]/warp\.write = \[1 1\]/' Huang_et_al_2013/start_seg.m \
    && grep -q 'warp.write = \[1 1\]' Huang_et_al_2013/start_seg.m \
    && patch -p1 -d Huang_et_al_2013 < /tmp/patches/mysegment-component-count.patch

# NIfTI and ANALYZE tools (Jimmy Shen, File Exchange 8797)
RUN wget -q "https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/8797/versions/28/download/zip" \
        -O /tmp/nifti_tools.zip \
    && unzip -q /tmp/nifti_tools.zip -d ./NIfTI_tools \
    && rm -f /tmp/nifti_tools.zip \
    && test -f ./NIfTI_tools/load_untouch_nii.m

# HArtMuT artefact source models + individual warp templates. Only the large
# model is stored in git-lfs; the pipeline uses HArtMuT_NYhead_small.mat.
RUN GIT_LFS_SKIP_SMUDGE=1 git clone --depth 1 https://github.com/harmening/hartmut.git HArtMuT

# weighted-SPHARM (Chung et al.)
RUN mkdir -p ./weighted-SPHARM \
    && wget -q https://pages.stat.wisc.edu/~mchung/research/amygdala/SPHARMsmooth2.m \
            -O ./weighted-SPHARM/SPHARMsmooth2.m \
    && wget -q https://pages.stat.wisc.edu/~mchung/softwares/weighted-SPHARM/SPHARMvectorize.m \
            -O ./weighted-SPHARM/SPHARMvectorize.m

# ATRA (ACPC alignment) - optional, license requires a manual download from
# https://www.nitrc.org/frs/download.php/10393/atra1.0_LinuxCentOS6.7.tar.gz
# Place atra1.0_LinuxCentOS6.7.tar.gz in the build context to include it.
RUN --mount=type=bind,source=.,target=/ctx \
    mkdir -p ./art \
    && if ls /ctx/atra1.0_*.tar.gz >/dev/null 2>&1; then \
         tar -xzf /ctx/atra1.0_*.tar.gz -C ./art \
         && echo "ATRA installed"; \
       else \
         echo "ATRA archive not found in build context, preprocessing (ACPC) will be unavailable"; \
       fi
ENV ARTHOME=${MRISEG_HOME}/art
ENV PATH=${ARTHOME}/bin:${PATH}

# FieldTrip compatibility folders for old MATLAB releases shadow builtins and
# break SPM12 (see README troubleshooting)
RUN find ./fieldtrip/compat -mindepth 1 -maxdepth 1 -type d -name "matlablt*" -exec rm -rf {} + \
    && find ./fieldtrip/external/spm12/external/fieldtrip/compat -mindepth 1 -maxdepth 1 \
            -type d -name "matlablt*" -exec rm -rf {} + 2>/dev/null || true

# ============================================================================
# Pipeline code from this repository (last: changes most often)
# ============================================================================
COPY start_segmentation.m start_segmentation_cedalion.m ./
COPY MRIsegmentation ./MRIsegmentation
COPY Schaefer2018_Parcellations ./Schaefer2018_Parcellations
COPY fieldtrip_mni_templates ./fieldtrip_mni_templates
COPY docker/mriseg_init.m ./docker/mriseg_init.m
COPY docker/mriseg /usr/local/bin/mriseg
RUN chmod 0755 /usr/local/bin/mriseg \
    && chown -R matlab:matlab ${MRISEG_HOME}

# mriseg_init.m (on MATLABPATH) sets up all paths; the wrapper calls it before
# every command. MATLABPATH keeps the entry scripts reachable for plain
# `matlab -batch` calls too.
ENV MATLABPATH=${MRISEG_HOME}:${MRISEG_HOME}/docker:${MRISEG_HOME}/MRIsegmentation

# MathWorks' reference Dockerfile force-enables usage telemetry (DDUX) here;
# it crashed MATLAB with heap corruption in the container, so it stays off.
ENV MW_DDUX_FORCE_ENABLE=false

USER matlab
WORKDIR ${MRISEG_HOME}

ENTRYPOINT ["tini", "--"]
CMD ["mriseg", "--help"]
