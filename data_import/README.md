# Data Import Tooling

The `data_import/import_dataset.py` helper automates the manual workflow of exec’ing into
the running PostgreSQL container and issuing `psql` commands. With a YAML file per
dataset, operators can recreate database refreshes without typing the individual steps.

## Prerequisites

- A PostgreSQL container is running (e.g. `descriptor-libraries-descriptor_database-1`) built
  from this repo’s `database/Dockerfile`, which installs the RDKit extension. Build it with
  `docker build -t descriptor-db -f database/Dockerfile .` and base your runtime container on
  that image so the extension and its dependencies are present.
- The CSV data and any schema SQL live inside that container. In production we bind a
  host directory to `/scratch` (for example `/PATH/TO/SHARED_SCRATCH:/scratch`) so the
  importer can read large CSV files and schema dumps from that mount. The demo stack uses
  `/docker-entrypoint-initdb.d/` for its bundled data.
- Python 3 and PyYAML are available on the host invoking the script.

## Usage

1. Copy one of the YAML examples in `data_import/configs/` and edit the paths/columns for
   your dataset.
2. Run the importer from the host:

   ```bash
   python data_import/import_dataset.py --config data_import/configs/<your-config>.yaml
   ```

The helper will drop/create the database (if requested), ensure extensions exist, build
or update property tables, `\copy` each CSV, and run any follow-up SQL you specify.

## YAML Configuration Overview

```yaml
container_name: descriptor-libraries-descriptor_database-1
database:
  name: library_db
  drop_database: true

extensions:
  - cube
  - rdkit

property_sets:
  dft:
    table_name: public.dft_data
    columns:
      - name: molecule_id
        type: text
        nullable: false
      - name: property
        type: text
        nullable: false
      - name: max
        type: numeric
      - name: min
        type: numeric
      - name: low_e
        type: numeric
      - name: boltz
        type: numeric
        csv: Boltz
    primary_key: [molecule_id, property]

schema_files:
  - /scratch/schema_dump.sql

tables:
  - name: public.molecule
    csv_path: /scratch/molecule.csv
    columns: [molecule_id, smiles, canonical_smiles, pca, umap]
    post_sql:
      - "UPDATE public.molecule SET mol = mol_from_smiles(smiles) WHERE mol IS NULL;"
      - "UPDATE public.molecule SET morganbv = morganbv_fp(mol) WHERE mol IS NOT NULL AND morganbv IS NULL;"

  - property_set: dft
    csv_path: /scratch/dft_properties.csv

post_import_sql:
  - "ANALYZE;"
```

### Key Sections

- **container_name**: Name of the running PostgreSQL container; all commands use
  `docker exec` against it.
- **database**: Target database options.
  - `name`: Database to populate.
  - `drop_database` (default `false`): Drop/recreate before loading.
  - `maintenance_db` / `owner` (optional): Override admin database or final owner.
- **extensions**: Extensions to ensure are installed (usually `cube` and `rdkit`).
- **schema_files**: SQL files (inside the container) executed with `psql -f` before imports.
- **property_sets**: Reusable definitions for property tables (`dft`, `ml`, `qikprop`, etc.).
  - `table_name`: Destination table. Defaults to `<property_name>_data` if omitted.
  - `columns`: Ordered column definitions (`name`, `type`, optional `nullable`, `default`,
    `csv`). Column order must match the CSV column order (the `csv` field is a reminder if
    the header label differs from the column name).
  - `primary_key`: Optional list of column names when you need a composite key.
  - `molecule_fk` (default `true` when `molecule_id` exists): Adds a foreign key back to
    `public.molecule(molecule_id)`; set `molecule_fk_on_delete: CASCADE` to cascade deletes.
  - `indexes` / `post_create_sql`: Additional SQL statements. Use
    `CREATE INDEX IF NOT EXISTS ...` to remain idempotent; `{table}` is replaced with the
    fully qualified table name.
- **tables**: Ordered list of CSV imports.
  - Use `name` + `columns` for explicit tables like `public.molecule` or `public.conformer`.
  - Use `property_set: <key>` to instantiate one of the reusable property definitions.
  - `pre_sql` / `post_sql`: Optional statements executed around that specific import.
  - `copy_options`: Custom `WITH (...)` clause for `\copy` (defaults to `FORMAT csv, HEADER true`).
- **post_import_sql**: Statements executed after all tables finish loading.

## Included Configurations

- `demo.yaml`: Demonstrates ingesting the repo’s demo CSVs plus a reusable DFT property set.
- `production.template.yaml`: Starting point for a production deployment that mirrors the
  descriptor libraries server. Fill in the database name, update CSV paths on shared
  storage (e.g., `/mnt/largestore1/...`), and expand property sets as needed.

## Legacy Notes

If you prefer to run everything by hand, keep a copy of your step-by-step notes alongside
these configs. The importer simply captures those same commands in an executable form.
