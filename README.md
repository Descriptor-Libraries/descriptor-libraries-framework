# Molecular Library Web Application

This repository provides a web application for displaying molecule library information. Originally developed for "kraken - Kolossal viRtual dAtabase for moleKular dEscriptors of orgaNophosphorus ligands", the platform has evolved into a flexible, multi-tenant system (one codebase serving multiple independent libraries) that can host multiple molecular databases with custom branding. You can see all available libraries at [descriptor-libraries.molssi.org](https://descriptor-libraries.molssi.org/).

The app is built with FastAPI (a Python web framework) for the backend and React (a JavaScript UI library) for the frontend, deployed as Docker containers that connect to a PostgreSQL database with the RDKit extension for chemistry operations, containing molecule and conformer information.

## Tech Stack

- **Backend**: FastAPI on Uvicorn (high-performance Python web server), SQLAlchemy (automap; reflects tables at runtime), Pydantic Settings (configuration management), Alembic (database migration tool)
- **Database**: PostgreSQL with the RDKit extension (adds chemistry-specific functions like structure searching)
- **Frontend**: React + Vite (fast build tool) with MUI (Material-UI components), Plotly (interactive graphs), Ketcher (molecule editor), built to static files and served by a small Node server (`server.js`)
- **Reverse proxy**: Traefik v2 for HTTPS certificates and path-based routing (directing URLs to the right service)
- **Depictions**: CDK Depict service for converting SMILES notation (text representation of molecules) to images
- **Containers**: Docker for packaging all services; Docker Compose for local/dev orchestration

## Architecture Overview

The application supports multiple molecular libraries through a namespace-based architecture:
- Each library gets its own URL path (e.g., `/acids`, `/bases`, `/kraken`)
- Custom branding and content per library (logos, color schemes, text)
- Shared backend infrastructure with library-specific API endpoints (e.g., `/api/acids`)
- Path-based routing through a reverse proxy (incoming requests are directed based on the URL path)

## Configuration & Paths

- UI path: set `VITE_BASE_URL=/<section>` (e.g., `/acids`).
- API path: set `API_PREFIX=<section>` so the backend serves `/api/<section>` (docs at `/api/<section>/docs`).

Running the Demo Application
=============================

A minimal demo version is available with 10 sample cyanoarene molecules to test the application locally.

## Quick Start

1. **Start the demo**:
   ```bash
   docker-compose up --build
   ```

2. **Access the application**:
   - **Web App**: http://localhost:9080/demo
   - **API Docs**: http://localhost:9080/api/demo/docs
   - **Traefik Dashboard**: http://localhost:9081
   - **Database**: localhost:5433 (user: postgres, password: postgres, db: demo_data)

3. **Stop the application**:
   ```bash
   docker-compose down
   ```

## What's Included in the Demo

- **10 Cyanoarene Molecules**: Representative sample with SMILES, 3D coordinates, PCA/UMAP projections
- **DFT Properties**: HOMO, LUMO, dipole moment, chemical hardness (η)
- **Chemical Search**: RDKit-powered substructure and similarity search
- **Interactive UI**: Browse, search, and visualize molecular data

The database automatically initializes with demo data on first run.

Running the Application for Production
======================================

Coming soon - Production deployment configurations are being updated.

Testing
=======

Coming soon - Test suite updates for the namespace-aware architecture are in progress.

Acknowledgements
----------------

This project is a collaboration between [The Molecular Sciences Software Institute](https://molssi.org/) and the [Center for Computer Assisted Synthesis](https://ccas.nd.edu/).

MolSSI is funded by the National Science Foundation OAC-1547580 and CHE-2136142.
C-CAS is funded by the National Science Foundation CHE–2202693.
