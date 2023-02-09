---- Set column pca in table pca to an array of all 4 pca component columns ----
ALTER TABLE pca ADD COLUMN pca float[];
UPDATE pca
SET pca = ARRAY[pca1, pca2, pca3, pca4];

---- Query to add molecule_id to umap table ----

ALTER TABLE new_data
RENAME TO umap;

ALTER TABLE umap ADD COLUMN molecule_id INT NOT NULL DEFAULT 0;

UPDATE umap t1
SET molecule_id = t2.molecule_id
FROM molecule t2
WHERE t1.smiles = t2.smiles;

-- check to make sure everything was transfered over
select molecule_id, smiles from umap
WHERE molecule_id IS NULL;

---- Query to add molecule_id to pca table ----

ALTER TABLE pca ADD COLUMN molecule_id INT NOT NULL DEFAULT 0;

UPDATE pca t1
SET molecule_id = t2.molecule_id
FROM molecule t2
WHERE t1.smiles = t2.smiles;

-- check to make sure everything was transfered over
select molecule_id, smiles from pca
WHERE molecule_id IS NULL;

---- Query to drop unnescesary columns from umap and pca tables ----

ALTER TABLE umap
DROP COLUMN umap1, DROP COLUMN umap2;

ALTER TABLE pca
DROP COLUMN pca1, DROP COLUMN pca2, DROP COLUMN pca3, DROP COLUMN pca4;