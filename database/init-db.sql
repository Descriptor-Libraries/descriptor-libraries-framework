-- Create database for demo
CREATE DATABASE demo_data;

-- Connect to demo database
\c demo_data;

-- Create extensions
CREATE EXTENSION IF NOT EXISTS cube WITH SCHEMA public;
CREATE EXTENSION IF NOT EXISTS rdkit WITH SCHEMA public;

-- Create molecule table
CREATE TABLE public.molecule (
    molecule_id text NOT NULL,
    smiles text,
    canonical_smiles text,
    mol public.mol,
    pca public.cube,
    umap public.cube,
    molecular_weight double precision,
    pat text,
    morganbv public.bfp
);

-- Create conformer table  
CREATE TABLE public.conformer (
    molecule_id text,
    coords double precision[],
    elements character varying[],
    basis text,
    method text,
    program text,
    cluster boolean,
    conformer_id SERIAL PRIMARY KEY
);

-- Create dft_data table
CREATE TABLE public.dft_data (
    molecule_id text,
    property text,
    boltz numeric,
    high_E numeric,
    low_E numeric
);

-- Set ownership
ALTER TABLE public.molecule OWNER TO postgres;
ALTER TABLE public.conformer OWNER TO postgres;
ALTER TABLE public.dft_data OWNER TO postgres;

-- Add constraints
ALTER TABLE ONLY public.molecule ADD CONSTRAINT molecule_pkey PRIMARY KEY (molecule_id);
ALTER TABLE ONLY public.molecule ADD CONSTRAINT molecule_smiles_key UNIQUE (smiles);

-- Add foreign key constraints
ALTER TABLE ONLY public.conformer ADD CONSTRAINT conformer_molecule_id_fkey FOREIGN KEY (molecule_id) REFERENCES public.molecule(molecule_id);
ALTER TABLE ONLY public.dft_data ADD CONSTRAINT dft_data_molecule_id_fkey FOREIGN KEY (molecule_id) REFERENCES public.molecule(molecule_id);

-- Create indexes
CREATE INDEX idx_molecule_molecule_id ON public.molecule USING btree (molecule_id);
CREATE INDEX idx_molecule_smiles ON public.molecule USING btree (smiles);
CREATE INDEX idx_molecule_pca ON public.molecule USING btree (pca);
CREATE INDEX idx_molecule_pca_knn ON public.molecule USING gist (pca);
CREATE INDEX idx_molecule_umap ON public.molecule USING btree (umap);
CREATE INDEX idx_molecule_umap_knn ON public.molecule USING gist (umap);
CREATE INDEX molecule_mol ON public.molecule USING gist (mol);
CREATE INDEX morganbv_idx ON public.molecule USING gist (morganbv);
CREATE INDEX mw_idx ON public.molecule USING btree (molecular_weight);