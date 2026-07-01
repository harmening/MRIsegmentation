# MRIsegmentation Dockerfile
#
# Build:   podman build -t mrisegmentation .
#
# Run segmentation:
#   podman run -it -v /path/to/data:/data mrisegmentation -batch "start_segmentation_cedalion('/data/T1.nii')"
# Interactive MATLAB:
#   podman run -it --shm-size=512M -p 8888:8888 -v /path/to/data:/data mrisegmentation -browser
#
# License: Requires MATLAB license configured for cloud use.
#          Individual and Campus-Wide licenses are already configured for cloud use.
#          See: https://www.mathworks.com/help/cloudcenter/ug/matlab-container-on-docker-hub.html

# ============================================================================
# Build arguments
# ============================================================================
ARG MATLAB_IMAGE=matlab:r2024b
ARG MATLAB_RELEASE=R2024b
ARG MATLAB_TOOLBOXES="Image_Processing_Toolbox Optimization_Toolbox Parallel_Computing_Toolbox Signal_Processing_Toolbox Statistics_and_Machine_Learning_Toolbox Wavelet_Toolbox"

# ============================================================================
# Base image
# ============================================================================
FROM mathworks/${MATLAB_IMAGE}

ARG MATLAB_RELEASE
ARG MATLAB_TOOLBOXES

LABEL maintainer="Nils Harmening"
LABEL description="MRI Segmentation Pipeline with SPM12, CAT12, FieldTrip, iso2mesh"

# Avoid interactive prompts
ENV DEBIAN_FRONTEND=noninteractive

# Switch to root for installation
USER root

ENV MATLAB_INSTALL_LOCATION=/opt/matlab/${MATLAB_RELEASE}
ENV PATH=${MATLAB_INSTALL_LOCATION}/bin:${PATH}

# ============================================================================
# Install system dependencies
# ============================================================================
RUN sudo apt-get update && sudo apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    curl \
    git \
    sudo \
    unzip \
    wget \
    xvfb \
    && rm -rf /var/lib/apt/lists/*


# ============================================================================
# Create working directory
# ============================================================================
ENV MRISEG_HOME=/opt/MRIsegmentation
WORKDIR ${MRISEG_HOME}

# ============================================================================
# Clone MRIsegmentation repository (or copy local files)
# ============================================================================
RUN git clone --depth 1 https://github.com/harmening/MRIsegmentation.git .

# ============================================================================
# Install FieldTrip
# ============================================================================
RUN git clone --depth 1 https://github.com/fieldtrip/fieldtrip.git

# ============================================================================
# Install SPM12 
# ============================================================================
#RUN wget -q https://www.fil.ion.ucl.ac.uk/spm/download/restricted/eldorado/spm12.zip && \

RUN --mount=type=bind,source=spm12.zip,target=/tmp/spm12.zip \
    rm -rf ./fieldtrip/external/spm12 && \
    unzip -q /tmp/spm12.zip -d ./fieldtrip/external/ 

# ============================================================================
# Build SPM12 MEX files
# ============================================================================
RUN cd ./fieldtrip/external/spm12/src && \
    make distclean && \
    make && make install && \
    make external-distclean && \
    make external && make external-install && \
    cd ../../../../

# ============================================================================
# Install CAT12 (version r2159 for compatibility with SPM12)
# ============================================================================
#RUN wget -q http://www.neuro.uni-jena.de/cat12/cat12_r2159.zip && \

RUN --mount=type=bind,source=cat12_r2159.zip,target=/tmp/cat12_r2159.zip \
    rm -rf ./fieldtrip/external/spm12/toolbox/cat12 && \
    unzip -q /tmp/cat12_r2159.zip -d ./fieldtrip/external/spm12/toolbox/ 

# ============================================================================
# Install iso2mesh
# ============================================================================
RUN git clone --depth 1 https://github.com/fangq/iso2mesh.git

# ============================================================================
# Install Andy's tools (Huang et al. 2013)
# ============================================================================
RUN --mount=type=bind,source=Huang_et_al_2013.zip,target=/tmp/Huang_et_al_2013.zip \
    unzip -q /tmp/Huang_et_al_2013.zip -d ./Huang_et_al_2013

# Apply patch: Make new_segment store the nonlinear warp
RUN sed -i 's/warp\.write = \[0 0\]/warp\.write = \[1 1\]/' Huang_et_al_2013/start_seg.m

# ============================================================================
# Install NIfTI tools
# ============================================================================
RUN wget -q "https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/8797/versions/28/download/zip" \
    -O nifti_tools.zip && \
    unzip -q nifti_tools.zip -d ./NIfTI_tools && \
    rm nifti_tools.zip || echo "NIfTI tools may need manual installation"

# ============================================================================
# Install weighted-SPHARM (optional)
# ============================================================================
RUN mkdir -p ./weighted-SPHARM && \
    wget -q https://pages.stat.wisc.edu/~mchung/research/amygdala/SPHARMsmooth2.m \
         -O ./weighted-SPHARM/SPHARMsmooth2.m && \
    wget -q https://pages.stat.wisc.edu/~mchung/softwares/weighted-SPHARM/SPHARMvectorize.m \
         -O ./weighted-SPHARM/SPHARMvectorize.m || \
    echo "weighted-SPHARM download failed - optional component"

# ============================================================================
# ATRA (Automatic Temporal Registration Algorithm) - OPTIONAL
# Due to license restrictions, must be downloaded manually from:
# https://www.nitrc.org/frs/download.php/10393/atra1.0_LinuxCentOS6.7.tar.gz
#
# To include ATRA:
# 1. Download and place atra1.0_LinuxCentOS6.7.tar.gz in build context
# 2. Uncomment the following lines:
# ============================================================================
RUN mkdir -p ./art
# COPY atra1.0_LinuxCentOS6.7.tar.gz ./art/
# RUN cd ./art && \
#     gunzip atra1.0_LinuxCentOS6.7.tar.gz && \
#     tar -xvf atra1.0_LinuxCentOS6.7.tar && \
#     rm -f atra1.0_*.tar

ENV ARTHOME=${MRISEG_HOME}/art
ENV PATH=${ARTHOME}/bin:${PATH}

# ============================================================================
# Install MATLAB using mpm (MathWorks Package Manager)
# https://github.com/mathworks-ref-arch/matlab-dockerfile
# Placed after other downloads so toolbox changes don't invalidate cache
# ============================================================================

RUN wget -q https://www.mathworks.com/mpm/glnxa64/mpm && \
    chmod +x mpm && \
    ./mpm install \
        --release=${MATLAB_RELEASE} \
        --destination=${MATLAB_INSTALL_LOCATION} \
        --products MATLAB ${MATLAB_TOOLBOXES} \
    || (cat /tmp/mathworks_root.log && false) && \
    rm -f mpm /tmp/mathworks_root.log 


# ============================================================================
# Fix FieldTrip compatibility issues (see README troubleshooting)
# ============================================================================
RUN find ./fieldtrip/compat -mindepth 1 -maxdepth 1 -type d -name "matlablt*" \
    -exec rm -rf {} + 2>/dev/null || true && \
    find ./fieldtrip/external/spm12/external/fieldtrip/compat -mindepth 1 -maxdepth 1 \
    -type d -name "matlablt*" -exec rm -rf {} + 2>/dev/null || true

# ============================================================================
# Copy startup.m for MATLAB path configuration
# ============================================================================
RUN mkdir -p /home/matlab/Documents/MATLAB
COPY docker/startup.m /home/matlab/Documents/MATLAB/startup.m

# Set MATLABPATH so paths are available in batch mode
# Note: startup.m is still needed for ft_defaults and other initialization
ENV MATLABPATH=${MRISEG_HOME}:${MRISEG_HOME}/MRIsegmentation:${MRISEG_HOME}/fieldtrip:${MRISEG_HOME}/iso2mesh:${MRISEG_HOME}/Huang_et_al_2013:${MRISEG_HOME}/NIfTI_tools:${MRISEG_HOME}/weighted-SPHARM:${MRISEG_HOME}/Schaefer2018_Parcellations

# ============================================================================
# Set permissions
# ============================================================================
RUN chown -R matlab:matlab ${MRISEG_HOME} && \
    chown -R matlab:matlab /home/matlab

# ============================================================================
# Switch to matlab user
# ============================================================================
USER matlab

# ============================================================================
# Data volume mount point
# ============================================================================
VOLUME ["/data"]

# ============================================================================
# Working directory for MATLAB
# ============================================================================
WORKDIR ${MRISEG_HOME}

# ============================================================================
# Use the MathWorks run.sh for license validation
# ============================================================================
ENTRYPOINT ["/bin/run.sh"]
CMD ["-help"]
