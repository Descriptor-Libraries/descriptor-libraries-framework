#!/bin/bash

# Initialize the shell for micromamba
eval "$(micromamba shell hook --shell bash)"

# Activate the environment. If you're using a different environment name other than "base", replace it.
micromamba activate base

# Navigate to app directory and start FastAPI using Uvicorn.
cd app
uvicorn main:app --reload --port 8080 --host 0.0.0.0
