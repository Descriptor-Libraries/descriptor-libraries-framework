# Molecular Library Web Application

This repository provides a web application for displaying molecule library information. Originally developed for "kraken - Kolossal viRtual dAtabase for moleKular dEscriptors of orgaNophosphorus ligands", the platform has evolved into a flexible, multi-tenant system (one codebase serving multiple independent libraries) that can host multiple molecular databases with custom branding. You can see all available libraries at [descriptor-libraries.molssi.org](https://descriptor-libraries.molssi.org/).

The app is built with FastAPI (a Python web framework) for the backend and React (a JavaScript UI library) for the frontend, deployed as Docker containers that connect to a PostgreSQL database with the RDKit extension for chemistry operations, containing molecule and conformer information.

## Quick Start

A minimal demo version is available with 10 sample cyanoarene molecules to test the application locally. To start, clone the repository

1. **Clone the repository**
   ```baseh
   git clone https://github.com/Descriptor-Libraries/descriptor-libraries-framework.git
   cd descriptor-libraries-framework
   ```
2. **Start the demo**:
   ```bash
   docker compose up
   ```

3. **Access the application**:
   - **Web App**: http://localhost/demo
   - **API Docs**: http://localhost/api/demo/docs
   - **Traefik Dashboard**: http://localhost:8080
   - **Database**: localhost:5433 (user: postgres, password: postgres, db: demo_data)

4. **Stop the application**:
   ```bash
   docker compose down
   ```

## What's Included in the Demo

- **10 Cyanoarene Molecules**: Representative sample with SMILES, 3D coordinates, PCA/UMAP projections
- **DFT Properties**: HOMO, LUMO, dipole moment, chemical hardness (η)
- **Chemical Search**: RDKit-powered substructure and similarity search
- **Interactive UI**: Browse, search, and visualize molecular data

The database automatically initializes with demo data on first run (the Postgres container creates `demo_data` and runs bundled init scripts to create the schema and load CSVs


## Tech Stack and Architecture

- **Backend**: FastAPI on Uvicorn (high-performance Python web server), SQLAlchemy, Pydantic Settings (configuration management)
- **Database**: PostgreSQL with the RDKit extension
- **Frontend**: React + Vite (build tool) with MUI (Material-UI components), Plotly (interactive graphs), Ketcher (molecule editor), built to static files and served by a small Node server (`server.js`)
- **Reverse proxy**: Traefik v2 for HTTPS certificates and path-based routing (directing URLs to the right service)
- **Depictions**: CDK Depict service for converting SMILES notation to images
- **Containers**: Docker for packaging all services; Docker Compose for local/dev orchestration

### How Namespaces Work
The application supports multiple molecular libraries through a **namespace-based architecture**. Each library is deployed as a separate instance with its own namespace:
- **URL Namespace**: Each library gets its own URL path (e.g., `/acids`, `/cyanoarenes`, `/demo`)
- **API Namespace**: Backend serves at `/api/<namespace>` (e.g., `/api/acids/molecules`)
- **Database Isolation**: Each library has its own database (e.g., `acids_data`, `demo_data`)
- **Custom Branding**: Library-specific logos, colors, and content

### Configuration Pattern
Libraries are configured through environment variables. The frontend uses `VITE_BASE_URL=/<namespace>` (like `/demo`), the backend uses `API_PREFIX=<namespace>` (like `demo`), and each gets its own database with `POSTGRES_DB=<namespace>_data` (like `demo_data`). This creates clean URL routing where frontends are served at `/<namespace>`, APIs at `/api/<namespace>`, and documentation at `/api/<namespace>/docs`.

### Service Components

The Docker Compose deployment consists of these services, each with specialized startup processes:

- **reverse-proxy (Traefik)**: Entrypoint starts Traefik with dashboard and routes requests to services based on URL
- **database (PostgreSQL + RDKit)**: Entrypoint initializes PostgreSQL; on first start creates the demo database and runs bundled SQL to set up tables and load demo data
- **backend (FastAPI)**: Entrypoint runs `prestart.sh`, which prepares the Python environment and starts the API with auto-reload
- **cdk-depict**: Entrypoint starts the depiction service used to render molecule images
- **frontend (React)**: Entrypoint runs `renameBase-dev.sh`, which updates the app's base path (e.g., `/demo`) and starts the dev server

### Frontend Base Path Handling

The `renameBase-dev.sh` script is essential for multi-tenant support. It replaces the placeholder `/base_url` in generated files with the value of `VITE_BASE_URL`, allowing the frontend to work under namespace paths like `/demo` instead of only at the site root `/`. After adjusting links and asset paths, it starts the development server so routes and static assets resolve correctly behind the reverse proxy.
Running the Demo Application
=============================




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
