#!/usr/bin/env bash
# Download the archives the Dockerfile bind-mounts from the build context.
# Run from the repository root:  bash docker/fetch-deps.sh
#
# ATRA (optional, ACPC alignment) has to be downloaded manually after accepting
# the license at
#   https://www.nitrc.org/frs/download.php/10393/atra1.0_LinuxCentOS6.7.tar.gz
# and placed next to this repository's Dockerfile.
set -euo pipefail
cd "$(dirname "$0")/.."

fetch() {
    local file="$1" url="$2"
    if [[ -s "$file" ]]; then
        echo "$file already present, skipping"
    else
        echo "downloading $file"
        wget -q --show-progress -O "$file" "$url"
    fi
}

fetch spm12.zip            https://www.fil.ion.ucl.ac.uk/spm/download/restricted/eldorado/spm12.zip
fetch cat12_r2159.zip      https://www.neuro.uni-jena.de/cat12/cat12_r2159.zip
fetch Huang_et_al_2013.zip https://www.parralab.org/segment/Huang_et_al_2013.zip
