FROM postgres:16.2

RUN apt update && \
    apt install -y --no-install-recommends postgresql-16-rdkit && \
    rm -rf /var/lib/apt/lists/*
