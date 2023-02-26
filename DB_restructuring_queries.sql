---- Add cube extension to DB ----
CREATE EXTENSION IF NOT EXISTS cube;

---- CREATE INDICES ON MOLECULE_ID ----
CREATE INDEX IF NOT EXISTS idx_molecule_molecule_id ON molecule (molecule_id);
CREATE INDEX IF NOT EXISTS idx_molecule_smiles ON molecule (smiles);
CREATE INDEX IF NOT EXISTS idx_molecule_umap_knn ON molecule USING gist (umap);
CREATE INDEX IF NOT EXISTS idx_molecule_pca ON molecule (pca);
CREATE INDEX IF NOT EXISTS idx_molecule_pca_knn ON molecule USING gist (pca);
CREATE INDEX IF NOT EXISTS idx_pca_smiles ON pca (smiles);

---- Set column pca in table molecule to a cube of all 4 pca component columns ----
ALTER TABLE molecule ADD COLUMN pca cube;
UPDATE molecule
SET pca = cube(ARRAY[pca.pca1, pca.pca2, pca.pca3, pca.pca4])
FROM pca
WHERE molecule.smiles = pca.smiles;

---- Query to make umap into cube data type on molecule table ----
ALTER TABLE molecule RENAME umap TO umap_point;
ALTER TABLE molecule ADD COLUMN umap cube;
UPDATE molecule
SET umap = cube(ARRAY[umap_point[0], umap_point[1]])
WHERE umap_point IS NOT NULL;

--- Drop umap_point ---
ALTER TABLE molecule DROP COLUMN umap_point;

---- Drop unneccessary tables ----
DROP TABLE new_data;
DROP TABLE pca;
DROP TABLE temporary;