#!/bin/sh

# Wait for DB to be ready
$python $(pwd)/app/backend_prestart.py

# Upgrade if needed
#alembic upgrade head

# Start FastAPI
cd app
uvicorn main:app --reload --port 8080 --host 0.0.0.0
