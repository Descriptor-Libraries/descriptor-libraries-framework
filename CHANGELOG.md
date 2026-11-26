# Changelog

Major changes in each release of the Descriptor Libraries framework: the shared
backend and frontend images that every library runs. Versions are
`vYEAR.MONTH.PATCH`.

## v2026.11.2 (2026-11-24)

### Fixed

- Molecule page: on a library without UMAP coordinates, choosing UMAP in the
  neighbors dropdown blanked the page. The UMAP option is now shown only when
  the library has UMAP neighbors.
- API: the neighbors endpoint returns an empty list when no molecule has
  coordinates of the requested type. It used to answer 500.

## v2026.11.1 (2026-11-16)

Replaces v2026.11.0, which was never deployed to production.

### Fixed

- RDKit pinned to 2024.03.5. The SMILES API matches RDKit's canonical SMILES
  against the `canonical_smiles` column stored in each library database. With
  RDKit 2025.03, tin and iron compounds could not be found by SMILES and were
  silently left out of batch exports. Upgrading RDKit later needs
  `canonical_smiles` recomputed in every library database at the same time.
- OpenBabel pinned to 3.1.1. Version 3.2 writes some ring bonds and valences
  differently in SDF conformer exports.
- The identifiers endpoint answered 500 on a library that stores PubChem IDs
  as text like `413605.0`. They are now read as numbers.

## v2026.11.0 (2026-11-13)

First shared release: one backend image and one frontend image for every
library.

### Added

- API: query a molecule by SMILES, and SMILES in batch downloads.
- API: `/version` endpoint reporting the version, git SHA and build date.
- Frontend: the footer shows the frontend and backend versions.
- Frontend: neighbor search components can be set per library with
  `neighbor_components` in `names.json`.
- Backend API tests, run with one Docker Compose command and by GitHub Actions
  on pushes and pull requests to `main`.

### Changed

- Neighbor search chooses its default components per library: the first three
  PCA components, or the first two if the library has only two.
- Backend dependencies are locked, and the backend starts in production mode.
- The frontend base image is pinned and installs from the lock file.
- Both images take `VERSION`, `GIT_SHA` and `BUILD_DATE` build arguments.
