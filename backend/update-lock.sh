#!/bin/bash
# Regenerate conda-linux-64.lock from environment.yaml.
#
# Run this after editing environment.yaml, or to pick up newer package builds
# within the pinned series. Commit environment.yaml and the lock file together,
# then rebuild the image and re-test before releasing.
set -euo pipefail
cd "$(dirname "$0")"

# Use the same micromamba image as the Dockerfile.
IMAGE=$(sed -n 's/^FROM //p' Dockerfile | head -1)

docker run --rm -v "$PWD/environment.yaml:/tmp/environment.yaml:ro" "$IMAGE" bash -c '
    micromamba install -y -q -n base -f /tmp/environment.yaml >/dev/null
    micromamba env export -n base --explicit --md5
' > conda-linux-64.lock.new

mv conda-linux-64.lock.new conda-linux-64.lock
echo "Wrote conda-linux-64.lock ($(grep -c '^https://' conda-linux-64.lock) packages)"
