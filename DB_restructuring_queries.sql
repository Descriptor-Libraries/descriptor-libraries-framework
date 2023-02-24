---- Add cube extension to DB ----
create extension cube;

---- Set column pca in table pca to a cube of all 4 pca component columns ----
ALTER TABLE pca ADD COLUMN pca cube;
UPDATE pca
SET pca = cube(ARRAY[pca1, pca2, pca3, pca4]);

---- Query to add molecule_id to umap table ----

ALTER TABLE new_data
RENAME TO umap;

ALTER TABLE umap ADD COLUMN molecule_id INT NOT NULL DEFAULT 0;

UPDATE umap t1
SET molecule_id = t2.molecule_id
FROM molecule t2
WHERE t1.smiles = t2.smiles;

-- check to make sure everything was transfered over
-- select molecule_id, smiles from umap
-- WHERE molecule_id IS NULL;

---- Query to add molecule_id to pca table ----

ALTER TABLE pca ADD COLUMN molecule_id INT NOT NULL DEFAULT 0;

UPDATE pca t1
SET molecule_id = t2.molecule_id
FROM molecule t2
WHERE t1.smiles = t2.smiles;

-- check to make sure everything was transfered over
-- select molecule_id, smiles from pca
-- WHERE molecule_id IS NULL;

---- Query to drop unnescesary columns from umap and pca tables ----

ALTER TABLE umap
DROP COLUMN umap1, DROP COLUMN umap2;

ALTER TABLE pca
DROP COLUMN pca1, DROP COLUMN pca2, DROP COLUMN pca3, DROP COLUMN pca4;

---- Creates distance function to calculate the distance between 2 arrays of floats ----
CREATE OR REPLACE FUNCTION distance(a float[], b float[])
RETURNS float AS $$
DECLARE
  n int := array_upper(a, 1);
  sum float := 0;
BEGIN
  FOR i IN 1..n LOOP
    sum := sum + (a[i] - b[i])^2;
  END LOOP;
  RETURN sqrt(sum);
END;
$$ LANGUAGE plpgsql;

---- Create new dft_data table and add only non-null dft_data into it ----

CREATE TABLE dft_data AS
SELECT molecule_id, smiles, dft_data
FROM molecule
WHERE dft_data IS NOT NULL;

-- check to make sure everything was transfered over
-- SELECT * FROM dft_data

---- Create new xtb_data table and add only non-null xtb_data into it ----

CREATE TABLE xtb_data AS
SELECT molecule_id, smiles, xtb_data
FROM molecule
WHERE xtb_data IS NOT NULL;

-- check to make sure everything was transfered over
-- SELECT * FROM xtb_data

---- Create new ml_data table and add only non-null ml_data into it ----

CREATE TABLE ml_data AS
SELECT molecule_id, smiles, ml_data
FROM molecule
WHERE ml_data IS NOT NULL;

-- check to make sure everything was transfered over
-- SELECT * FROM ml_data

---- Create new xtb_ni_data table and add only non-null xtb_ni_data into it ----

CREATE TABLE xtb_ni_data AS
SELECT molecule_id, smiles, xtb_ni_data
FROM molecule
WHERE xtb_ni_data IS NOT NULL;

-- check to make sure everything was transfered over
-- SELECT * FROM xtb_ni_data