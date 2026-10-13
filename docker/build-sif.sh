#!/usr/bin/env bash
# Convert the podman/docker image into an Apptainer/Singularity image file.
#
#   bash docker/build-sif.sh [image] [output.sif]
#
# Defaults: image localhost/mrisegmentation:latest, output ./mriseg.sif.
# Apptainer runs the container as you, with your home directory, hostname and
# network, which is what lets MATLAB inside reuse the MathWorks sign-in cached
# by the MATLAB installed on your host (online licensing, Individual licenses).
set -euo pipefail

image="${1:-localhost/mrisegmentation:latest}"
out="${2:-mriseg.sif}"
engine=podman
command -v podman >/dev/null 2>&1 || engine=docker
builder=apptainer
command -v apptainer >/dev/null 2>&1 || builder=singularity

tar="$(mktemp --tmpdir mriseg-XXXXXX.tar)"
trap 'rm -f "$tar"' EXIT

echo "exporting $image with $engine (about 12 GB) ..."
"$engine" save "$image" -o "$tar"
echo "building $out with $builder ..."
"$builder" build --force "$out" "docker-archive://$tar"
ls -lh "$out"
