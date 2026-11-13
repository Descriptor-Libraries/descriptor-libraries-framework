#!/bin/bash

# Initialize the shell for micromamba
eval "$(micromamba shell hook --shell bash)"

# Activate the environment. If you're using a different environment name other than "base", replace it.
micromamba activate base

# Navigate to app directory and start FastAPI using Uvicorn.
cd app

# Production start: no auto-reload. `exec` replaces this shell with uvicorn so
# that it receives stop signals directly and shuts down cleanly.
# - UVICORN_WORKERS sets the number of worker processes (default 1).
# - Extra arguments are passed on to uvicorn, e.g. `bash prestart.sh --reload`
#   for development.
exec uvicorn main:app --port 8080 --host 0.0.0.0 --workers "${UVICORN_WORKERS:-1}" "$@"
