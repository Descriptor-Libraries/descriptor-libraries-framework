# micromamba is pinned by version and digest so the base image cannot change underneath us.
FROM mambaorg/micromamba:2.10.0@sha256:e1476182cd40eb1f7c0cfae4adec9bfd78f3d7f05c61368306abfdc8939c1198

# Install exactly the packages recorded in the lock file (see environment.yaml
# and update-lock.sh for how it is produced).
COPY conda-linux-64.lock /tmp/conda-tmp/conda-linux-64.lock

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN micromamba install -y -n base -f /tmp/conda-tmp/conda-linux-64.lock && \
    micromamba clean --all --yes

WORKDIR /app/
ADD ./app /app/
# --no-deps: every dependency comes from the lock file, never from PyPI at build time.
RUN pip install --no-deps -e .
