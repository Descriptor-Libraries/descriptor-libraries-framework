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

### Library Settings in `names.json`

Each library's `public/brand/names.json` holds its name, tagline and title, plus these optional keys:

| Key | Meaning | When left out |
|---|---|---|
| `default_smiles` | Structure pre-filled on the search page | `C=C` |
| `default_display` | Projection on the home page plot: `"pca"` or `"umap"` | `umap` |
| `neighbor_components` | Principal components the neighbor search measures distance in, for example `"1,2"` | The backend uses the first three, or the first two if the library has only two |

Changing `neighbor_components` changes which molecules are listed as neighbors.

### Service Components

The Docker Compose deployment consists of these services, each with specialized startup processes:

- **reverse-proxy (Traefik)**: Entrypoint starts Traefik with dashboard and routes requests to services based on URL
- **database (PostgreSQL + RDKit)**: Entrypoint initializes PostgreSQL; on first start creates the demo database and runs bundled SQL to set up tables and load demo data
- **backend (FastAPI)**: Entrypoint runs `prestart.sh`, which prepares the Python environment and starts the API (the demo Compose file passes `--reload` for auto-reload; production runs it without)
- **cdk-depict**: Entrypoint starts the depiction service used to render molecule images
- **frontend (React)**: Entrypoint runs `renameBase-dev.sh`, which updates the app's base path (e.g., `/demo`) and starts the dev server

### Frontend Base Path Handling

The `renameBase-dev.sh` script is essential for multi-tenant support. It replaces the placeholder `/base_url` in generated files with the value of `VITE_BASE_URL`, allowing the frontend to work under namespace paths like `/demo` instead of only at the site root `/`. After adjusting links and asset paths, it starts the development server so routes and static assets resolve correctly behind the reverse proxy.
Running the Demo Application
=============================




Running the Application for Production
======================================

The production deployment mirrors the demo stack but keeps the PostgreSQL database
persistent and mounts a shared scratch directory for imports. Key differences:

- **Database image**: Build the Postgres + RDKit image that lives in `database/Dockerfile`
  and tag it (e.g., `docker build -t rdkit-postgres -f database/Dockerfile .`).
- **Persistent volumes**: Instead of the demo `tmpfs` database, mount the production
  storage so data survives restarts:

  ```yaml
  descriptor_database:
    image: rdkit-postgres
    shm_size: "4gb"
    volumes:
      - "/PATH/TO/PERSISTENT_DB:/var/lib/postgresql/data"
      - "/PATH/TO/SHARED_SCRATCH:/scratch"
  ```

  The first volume keeps the Postgres data directory on durable storage (e.g.,
  `/mnt/largestore1/descriptor-libraries/descriptor_postgres`). The `/scratch` mount makes
  large CSVs and schema dumps available inside the container for the import tooling (e.g.,
  `/mnt/largestore1/janash/ccas`). Adjust the host paths to match your environment.
- **Backups and imports**: Use the `data_import/import_dataset.py` helper (with the YAML
  configs under `data_import/configs/`) to load new libraries or refresh an existing one.

- **Backend image**: Build from a clean, tagged checkout and pass the version details, which
  the API then reports at `/api/<namespace>/version`:

  ```bash
  docker build -t descriptor-backend:v2026.11.0 \
    --build-arg VERSION=v2026.11.0 \
    --build-arg GIT_SHA="$(git rev-parse --short HEAD)" \
    --build-arg BUILD_DATE="$(date -u +%Y-%m-%dT%H:%M:%SZ)" \
    backend
  ```

  Without these arguments the endpoint reports `dev` / `unknown`.

Other services (reverse proxy, backend, frontend, depiction) run as in the demo compose
file. Ensure Traefik routes remain consistent with the namespaces you expose publicly.

Testing
=======

The backend API tests need only Docker. One command builds the backend, loads two
small test libraries, runs the tests and returns their exit code:

```bash
docker compose -f docker-compose.test.yml up --build --abort-on-container-exit --exit-code-from tests
docker compose -f docker-compose.test.yml down -v
```

No ports are published on the host, so this can run next to the demo or any other stack.

The two libraries are built from `tests/demo_data` and served by the same backend image:

| Library | Molecule IDs | `pubchem_id` column | Data |
|---|---|---|---|
| `demo` | text (`aryl-1` ...) | no | DFT |
| `numeric` (`tests/fixtures/numeric-ids.sql`) | integer | yes | DFT and ML |

To run the tests against an API that is already running, such as the demo above, set its
address (tests that need the second library are skipped):

```bash
pip install pytest requests
API_URL=http://localhost/api/demo pytest tests
```

Tests marked `xfail` describe requests the API currently answers with an error; see the
notes at the top of `tests/test_backend.py`.

The frontend has a lint check. It fails on errors and lists warnings:

```bash
cd frontend && yarn install --frozen-lockfile && npm run lint
```

Acknowledgements
----------------

This project is a collaboration between [The Molecular Sciences Software Institute](https://molssi.org/) and the [Center for Computer Assisted Synthesis](https://ccas.nd.edu/).

MolSSI is funded by the National Science Foundation OAC-1547580 and CHE-2136142.
C-CAS is funded by the National Science Foundation CHE–2202693.
