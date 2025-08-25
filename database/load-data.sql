-- Load demo data into tables
\c demo_data;

-- Copy data from CSV files
\copy public.molecule (molecule_id, smiles, canonical_smiles, pca, umap) FROM '/docker-entrypoint-initdb.d/csv/molecule.csv' WITH (FORMAT csv, HEADER);
\copy public.conformer (molecule_id, elements, coords, basis, program, method, cluster) FROM '/docker-entrypoint-initdb.d/csv/conformer.csv' WITH (FORMAT csv, HEADER);  
\copy public.dft_data (molecule_id, property, boltz, high_E, low_E) FROM '/docker-entrypoint-initdb.d/csv/dft_data.csv' WITH (FORMAT csv, HEADER);

-- Update RDKit molecule column
UPDATE public.molecule SET mol = mol_from_smiles(smiles::cstring);

-- Update morgan fingerprints
UPDATE public.molecule SET morganbv = morganbv_fp(mol);

-- Calculate molecular weights (if not provided)
UPDATE public.molecule SET molecular_weight = mol_amw(mol) WHERE molecular_weight IS NULL;