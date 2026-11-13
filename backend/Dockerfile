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

# What this image contains. Pass these at build time; they are reported by
# GET /api/<library>/version and recorded as image labels. They come last so
# that changing them does not invalidate the dependency layer above.
ARG VERSION=dev
ARG GIT_SHA=unknown
ARG BUILD_DATE=unknown
ENV APP_VERSION=${VERSION} \
    APP_GIT_SHA=${GIT_SHA} \
    APP_BUILD_DATE=${BUILD_DATE}
LABEL org.opencontainers.image.version=${VERSION} \
      org.opencontainers.image.revision=${GIT_SHA} \
      org.opencontainers.image.created=${BUILD_DATE} \
      org.opencontainers.image.source="https://github.com/Descriptor-Libraries/descriptor-libraries-framework" \
      org.opencontainers.image.url="https://github.com/Descriptor-Libraries/descriptor-libraries-framework" \
      org.opencontainers.image.title="descriptor-backend" \
      org.opencontainers.image.licenses="BSD-3-Clause" \
      org.opencontainers.image.description="Descriptor Libraries backend API (FastAPI + RDKit)"
