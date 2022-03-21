#!/bin/sh

# Wait for DB to be ready
python /app/app/backend_prestart.py

# Upgrade if needed
alembic upgrade head
