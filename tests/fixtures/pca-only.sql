-- Third test library: a copy of the numeric library without UMAP coordinates,
-- the shape of a library that only has PCA.
CREATE DATABASE test_pca_only TEMPLATE test_numeric;
\c test_pca_only;

UPDATE public.molecule SET umap = NULL;
