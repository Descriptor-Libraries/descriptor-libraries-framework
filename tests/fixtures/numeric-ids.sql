-- Second test library, built from the same demo CSV files, with the shape most
-- production libraries have: integer molecule IDs, a pubchem_id column,
-- ML data next to the DFT data and more than two PCA dimensions. The demo library (string IDs, no pubchem_id)
-- covers the other shape.
CREATE DATABASE test_numeric;
\c test_numeric;

CREATE EXTENSION IF NOT EXISTS cube WITH SCHEMA public;
CREATE EXTENSION IF NOT EXISTS rdkit WITH SCHEMA public;

CREATE TABLE public.molecule (
    molecule_id integer PRIMARY KEY,
    smiles text UNIQUE,
    canonical_smiles text,
    mol public.mol,
    pca public.cube,
    umap public.cube,
    molecular_weight double precision,
    pat text,
    morganbv public.bfp,
    pubchem_id integer
);

CREATE TABLE public.conformer (
    molecule_id integer REFERENCES public.molecule(molecule_id),
    coords double precision[],
    elements character varying[],
    basis text,
    method text,
    program text,
    cluster boolean,
    conformer_id SERIAL PRIMARY KEY
);

CREATE TABLE public.dft_data (
    molecule_id integer REFERENCES public.molecule(molecule_id),
    property text,
    boltz numeric,
    high_E numeric,
    low_E numeric
);

CREATE TABLE public.ml_data (
    molecule_id integer REFERENCES public.molecule(molecule_id),
    property text,
    boltz numeric,
    high_E numeric,
    low_E numeric
);

-- Load the CSV files as text first, then turn "aryl-7" into 7.
CREATE TEMP TABLE molecule_in (molecule_id text, smiles text, canonical_smiles text, pca public.cube, umap public.cube);
CREATE TEMP TABLE conformer_in (molecule_id text, elements character varying[], coords double precision[], basis text, program text, method text, cluster boolean);
CREATE TEMP TABLE dft_in (molecule_id text, property text, boltz numeric, high_E numeric, low_E numeric);

\copy molecule_in FROM '/docker-entrypoint-initdb.d/csv/molecule.csv' WITH (FORMAT csv, HEADER);
\copy conformer_in FROM '/docker-entrypoint-initdb.d/csv/conformer.csv' WITH (FORMAT csv, HEADER);
\copy dft_in FROM '/docker-entrypoint-initdb.d/csv/dft_data.csv' WITH (FORMAT csv, HEADER);

-- The demo data has 2 PCA dimensions. This library gets 4, so that requests for
-- three components can be tested: dimensions 3 and 4 are made up from the UMAP values.
INSERT INTO public.molecule (molecule_id, smiles, canonical_smiles, pca, umap)
    SELECT replace(molecule_id, 'aryl-', '')::integer, smiles, canonical_smiles,
           cube(ARRAY[pca->1, pca->2, (umap->1) / 3, (umap->2) / 3]), umap
    FROM molecule_in;

INSERT INTO public.conformer (molecule_id, elements, coords, basis, program, method, cluster)
    SELECT replace(molecule_id, 'aryl-', '')::integer, elements, coords, basis, program, method, cluster
    FROM conformer_in ORDER BY replace(molecule_id, 'aryl-', '')::integer;

INSERT INTO public.dft_data
    SELECT replace(molecule_id, 'aryl-', '')::integer, property, boltz, high_E, low_E FROM dft_in;

-- ML data for molecules 1 to 5 only, so "no data of this type" can be tested on the rest.
INSERT INTO public.ml_data
    SELECT replace(molecule_id, 'aryl-', '')::integer, property, boltz, high_E, low_E FROM dft_in
    WHERE replace(molecule_id, 'aryl-', '')::integer <= 5;

UPDATE public.molecule SET mol = mol_from_smiles(smiles::cstring);
UPDATE public.molecule SET morganbv = morganbv_fp(mol);
UPDATE public.molecule SET molecular_weight = mol_amw(mol);

-- A made-up CID for molecule 1, to prove the value is read from the database
-- and not looked up on PubChem. The other rows are left empty.
UPDATE public.molecule SET pubchem_id = 999000001 WHERE molecule_id = 1;

CREATE INDEX idx_molecule_pca_knn ON public.molecule USING gist (pca);
CREATE INDEX idx_molecule_umap_knn ON public.molecule USING gist (umap);
CREATE INDEX molecule_mol ON public.molecule USING gist (mol);
CREATE INDEX morganbv_idx ON public.molecule USING gist (morganbv);
