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