FROM mambaorg/micromamba

COPY environment.yaml* /tmp/conda-tmp/

ARG MAMBA_DOCKERFILE_ACTIVATE=1

RUN micromamba install -f /tmp/conda-tmp/environment.yaml && \
    micromamba clean --all --yes

WORKDIR /app/
ADD ./app /app/
RUN pip install -e .