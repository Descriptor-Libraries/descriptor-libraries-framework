"""
API tests for the backend.

They run against two small libraries built from tests/demo_data (see conftest.py
and docker-compose.test.yml):

    demo     string IDs ("aryl-1" ... "aryl-10"), no pubchem_id column, DFT data only
    numeric  integer IDs (1 ... 10), a pubchem_id column, ML data for molecules 1-5

Both have 10 cyanoarenes, one conformer each, and 2 PCA and 2 UMAP dimensions.

Tests marked `xfail(strict=True)` describe what the API should do where it
currently answers 500. They are reminders: when one of those is fixed the test
starts passing, pytest reports it, and the mark can be removed.
"""
import csv
import io

import pytest

# The same molecule written three ways: as stored, canonical, and a third valid form.
ARYL_4_STORED = "N#CC1=CC=NC=C1"
ARYL_4_CANONICAL = "N#Cc1ccncc1"
ARYL_4_OTHER = "c1cc(C#N)ccn1"
ARYL_4_INCHIKEY = "GPHQHTOMRSGBNZ-UHFFFAOYSA-N"

ARYL_1_CANONICAL = "Cc1cc(C#N)ccn1"
DFT_PROPERTIES = {"HOMO", "LUMO", "dipole", "η"}

known_500 = pytest.mark.xfail(strict=True, reason="currently answers 500")


@pytest.fixture(params=["demo", "numeric"])
def lib(request):
    """Run a test once per library. `lib.id(4)` is that library's ID for molecule 4."""
    api = request.getfixturevalue(request.param)
    api.name = request.param
    api.id = (lambda n: f"aryl-{n}") if request.param == "demo" else (lambda n: n)
    return api


# --- Service ---------------------------------------------------------------

def test_version(lib):
    r = lib.get("/version")
    assert r.status_code == 200
    assert set(r.json()) == {"version", "git_sha", "build_date"}


def test_v2_prefix_serves_the_same_api(lib):
    assert lib.get("/v2/molecules/first_id").json() == lib.get("/molecules/first_id").json()


# --- Molecules: string and integer IDs --------------------------------------

def test_first_id(lib):
    assert lib.get("/molecules/first_id").json() == {"first_molecule_id": lib.id(1)}


def test_first_id_keeps_its_type(demo, numeric):
    assert isinstance(demo.get("/molecules/first_id").json()["first_molecule_id"], str)
    assert isinstance(numeric.get("/molecules/first_id").json()["first_molecule_id"], int)


@pytest.mark.parametrize("n", [1, 4, 10])
def test_get_molecule(lib, n):
    r = lib.get(f"/molecules/{lib.id(n)}")
    assert r.status_code == 200
    body = r.json()
    assert body["molecule_id"] == lib.id(n)
    assert len(body["conformers_id"]) == 1
    assert body["molecular_weight"] > 100
    assert body["compound_name"] is None


def test_get_molecule_smiles(lib):
    assert lib.get(f"/molecules/{lib.id(4)}").json()["smiles"] == ARYL_4_STORED


@known_500
@pytest.mark.parametrize("missing", ["aryl-999", "999", "nope"])
def test_get_missing_molecule_is_404(lib, missing):
    assert lib.get(f"/molecules/{missing}").status_code == 404


# --- Data types ---------------------------------------------------------------

def test_data_types(demo, numeric):
    assert demo.get("/molecules/data_types").json() == {"available_types": ["dft_data"]}
    assert sorted(numeric.get("/molecules/data_types").json()["available_types"]) == ["dft_data", "ml_data"]


def test_molecule_data_types(demo, numeric):
    assert demo.get("/molecules/aryl-1/data_types").json() == {"dft_data": True}
    assert numeric.get("/molecules/1/data_types").json() == {"dft_data": True, "ml_data": True}
    assert numeric.get("/molecules/7/data_types").json() == {"dft_data": True, "ml_data": False}


# --- Data by ID ---------------------------------------------------------------

def test_data_by_id(lib):
    r = lib.get(f"/molecules/data/{lib.id(1)}", data_type="dft")
    assert r.status_code == 200
    rows = r.json()
    assert {row["property"] for row in rows} == DFT_PROPERTIES
    assert {row["molecule_id"] for row in rows} == {lib.id(1)}
    homo = next(row for row in rows if row["property"] == "HOMO")
    assert float(homo["boltz"]) == pytest.approx(-0.33537)


def test_data_by_id_defaults_to_ml(numeric):
    r = numeric.get("/molecules/data/1")
    assert r.status_code == 200
    assert {row["property"] for row in r.json()} == DFT_PROPERTIES


def test_data_by_id_without_that_data_is_empty(numeric):
    assert numeric.get("/molecules/data/7").json() == []


def test_data_by_id_invalid_type(lib):
    assert lib.get(f"/molecules/data/{lib.id(1)}", data_type="bad").status_code == 400


@known_500
def test_data_by_id_type_the_library_lacks(demo):
    # The demo library has no ml_data table at all.
    assert demo.get("/molecules/data/aryl-1", data_type="ml").status_code == 404


# --- Data by SMILES (the endpoints added on the smiles-query branch) ----------

@pytest.mark.parametrize("smiles", [ARYL_4_STORED, ARYL_4_CANONICAL, ARYL_4_OTHER])
def test_data_by_smiles_any_spelling(lib, smiles):
    r = lib.get("/molecules/data", smiles=smiles, data_type="dft")
    assert r.status_code == 200
    rows = r.json()
    assert {row["molecule_id"] for row in rows} == {lib.id(4)}
    assert {row["property"] for row in rows} == DFT_PROPERTIES


def test_data_by_smiles_matches_data_by_id(lib):
    by_smiles = lib.get("/molecules/data", smiles=ARYL_4_CANONICAL, data_type="dft").json()
    by_id = lib.get(f"/molecules/data/{lib.id(4)}", data_type="dft").json()
    assert by_smiles == by_id


def test_data_by_smiles_unknown_molecule(lib):
    r = lib.get("/molecules/data", smiles="CCCCCC", data_type="dft")
    assert r.status_code == 404
    assert r.json()["detail"] == "Molecule not found"


def test_data_by_smiles_invalid_smiles(lib):
    assert lib.get("/molecules/data", smiles="not_a_smiles", data_type="dft").status_code == 400


def test_data_by_smiles_type_missing_for_molecule(demo, numeric):
    # Default data type is ml: the demo library has none, numeric has none for molecule 7.
    assert demo.get("/molecules/data", smiles=ARYL_1_CANONICAL).status_code == 404
    assert numeric.get("/molecules/data", smiles="N#Cc1ccccn1").status_code == 404
    assert numeric.get("/molecules/data", smiles=ARYL_1_CANONICAL).status_code == 200


# --- Export -----------------------------------------------------------------------

def _csv_rows(response):
    return list(csv.DictReader(io.StringIO(response.text)))


def test_batch_export_by_ids_json(lib):
    ids = f"{lib.id(1)},{lib.id(2)}"
    r = lib.get("/molecules/data/export/batch", molecule_ids=ids, data_type="dft", return_type="json")
    assert r.status_code == 200
    rows = r.json()
    assert [row["molecule_id"] for row in rows] == [lib.id(1), lib.id(2)]
    assert rows[0]["HOMO_boltz"] == pytest.approx(-0.33537)


def test_batch_export_by_ids_csv(lib):
    ids = f"{lib.id(1)},{lib.id(2)}"
    r = lib.get("/molecules/data/export/batch", molecule_ids=ids, data_type="dft", context="substructure")
    assert r.status_code == 200
    assert r.headers["content-type"].startswith("text/csv")
    assert f"dft_{lib.id(1)}_2_substructure.csv" in r.headers["content-disposition"]
    rows = _csv_rows(r)
    assert [row["molecule_id"] for row in rows] == [str(lib.id(1)), str(lib.id(2))]
    assert "HOMO_boltz" in rows[0] and "smiles" in rows[0]


def test_batch_export_by_smiles(lib):
    smiles = f"{ARYL_4_CANONICAL}, N#Cc1ccccn1"
    r = lib.get("/molecules/data/export/batch", molecule_smiles=smiles, data_type="dft", return_type="json")
    assert r.status_code == 200
    assert sorted(str(row["molecule_id"]) for row in r.json()) == sorted([str(lib.id(4)), str(lib.id(7))])


def test_batch_export_by_smiles_matches_by_ids(lib):
    by_smiles = lib.get("/molecules/data/export/batch", molecule_smiles=ARYL_4_CANONICAL,
                        data_type="dft", return_type="json").json()
    by_ids = lib.get("/molecules/data/export/batch", molecule_ids=str(lib.id(4)),
                     data_type="dft", return_type="json").json()
    assert by_smiles == by_ids


@pytest.mark.parametrize(
    "params",
    [
        {},
        {"molecule_ids": "1", "molecule_smiles": "CC"},
        {"molecule_ids": "1", "data_type": "bad"},
        {"molecule_ids": "1", "return_type": "xml"},
        {"molecule_ids": "1", "context": "bad"},
    ],
)
def test_batch_export_invalid_parameters(lib, params):
    params.setdefault("data_type", "dft")
    assert lib.get("/molecules/data/export/batch", **params).status_code == 400


def test_single_export_csv(lib):
    r = lib.get(f"/molecules/data/export/{lib.id(1)}", data_type="dft")
    assert r.status_code == 200
    assert f"{lib.id(1)}_dft.csv" in r.headers["content-disposition"]
    rows = _csv_rows(r)
    assert len(rows) == 1 and rows[0]["molecule_id"] == str(lib.id(1))


def test_single_export_no_content(numeric):
    assert numeric.get("/molecules/data/export/7", data_type="ml").status_code == 204


# --- Substructure search ---------------------------------------------------------

def test_search_matches_all(lib):
    r = lib.get("/molecules/search/", substructure="C#N")
    assert r.status_code == 200
    assert {m["molecule_id"] for m in r.json()} == {lib.id(n) for n in range(1, 11)}


@pytest.mark.parametrize(
    "substructure, expected",
    [("c1ccsc1", {9}), ("c1ccoc1", {10}), ("B", {6}), ("CO", {2, 6, 8}), ("P", set())],
)
def test_search_hits(lib, substructure, expected):
    r = lib.get("/molecules/search/", substructure=substructure)
    assert r.status_code == 200
    assert {m["molecule_id"] for m in r.json()} == {lib.id(n) for n in expected}


def test_search_paging(lib):
    everything = [m["molecule_id"] for m in lib.get("/molecules/search/", substructure="C#N", limit=100).json()]
    page_1 = [m["molecule_id"] for m in lib.get("/molecules/search/", substructure="C#N", skip=0, limit=4).json()]
    page_2 = [m["molecule_id"] for m in lib.get("/molecules/search/", substructure="C#N", skip=4, limit=4).json()]
    assert len(page_1) == 4 and len(page_2) == 4
    assert not set(page_1) & set(page_2)
    assert set(page_1 + page_2) <= set(everything)


@pytest.mark.parametrize("substructure", ["C=XX", "C==C", "not a smiles"])
def test_search_invalid_smiles(lib, substructure):
    assert lib.get("/molecules/search/", substructure=substructure).status_code == 400


# --- Neighbors and dimensions -------------------------------------------------------

@pytest.mark.parametrize("type", ["pca", "umap", "PCA", "UMAP"])
def test_neighbors(lib, type):
    r = lib.get(f"/molecules/{lib.id(1)}/neighbors/", type=type, components="1,2", limit=5)
    assert r.status_code == 200
    rows = r.json()
    assert len(rows) == 5
    assert all(len(row["components"]) == 2 for row in rows)
    distances = [row["dist"] for row in rows]
    assert distances == sorted(distances)


def test_pca_neighbors_order(lib):
    rows = lib.get(f"/molecules/{lib.id(1)}/neighbors/", type="pca", components="1,2", limit=3).json()
    assert rows[0]["molecule_id"] == lib.id(2)
    assert rows[0]["dist"] == pytest.approx(0.2379220, abs=1e-6)


def test_neighbors_single_component(lib):
    rows = lib.get(f"/molecules/{lib.id(1)}/neighbors/", type="pca", components="1", limit=3).json()
    assert [len(row["components"]) for row in rows] == [1, 1, 1]


def test_neighbors_skip_zero_includes_the_molecule_itself(lib):
    # skip defaults to 1, which leaves the molecule itself out. The frontend
    # sends skip=0, so its first result is the molecule that was searched for.
    path = f"/molecules/{lib.id(1)}/neighbors/"
    default = lib.get(path, type="pca", components="1,2", limit=3).json()
    from_zero = lib.get(path, type="pca", components="1,2", skip=0, limit=4).json()
    assert from_zero[0]["molecule_id"] == lib.id(1)
    assert from_zero[0]["dist"] == 0
    assert from_zero[1:] == default


def test_neighbors_paging(lib):
    path = f"/molecules/{lib.id(1)}/neighbors/"
    all_six = [m["molecule_id"] for m in lib.get(path, type="pca", components="1,2", skip=0, limit=6).json()]
    second_page = [m["molecule_id"] for m in lib.get(path, type="pca", components="1,2", skip=3, limit=3).json()]
    assert second_page == all_six[3:]


@pytest.mark.parametrize(
    "params",
    [
        {"type": "bad", "components": "1"},
        {"type": "pca", "components": "1,2,3"},
        {"type": "pca", "components": "3"},
        {"type": "umap", "components": "3"},
        # With no components the backend asks for PCA 1,2,3; these libraries have 2.
        {"type": "pca"},
    ],
)
def test_neighbors_invalid_parameters(lib, params):
    assert lib.get(f"/molecules/{lib.id(1)}/neighbors/", **params).status_code == 400


def test_neighbors_default_umap_components(lib):
    assert lib.get(f"/molecules/{lib.id(1)}/neighbors/", type="umap", limit=2).status_code == 200


@pytest.mark.parametrize("type", ["pca", "umap"])
def test_dimensions(lib, type):
    r = lib.get("/molecules/dimensions/", type=type, components="1,2", limit=100)
    assert r.status_code == 200
    rows = r.json()
    assert {row["molecule_id"] for row in rows} == {lib.id(n) for n in range(1, 11)}
    assert all(len(row["components"]) == 2 and row["type"] == type for row in rows)


def test_dimensions_values(lib):
    rows = lib.get("/molecules/dimensions/", type="pca", components="1,2", limit=100).json()
    first = next(row for row in rows if row["molecule_id"] == lib.id(1))
    assert first["components"] == pytest.approx([-4.568498259542651, -1.382468555091367])


@pytest.mark.parametrize(
    "params", [{"type": "bad", "components": "1"}, {"type": "pca", "components": "3"}, {"type": "pca"}]
)
def test_dimensions_invalid_parameters(lib, params):
    assert lib.get("/molecules/dimensions/", **params).status_code == 400


@known_500
def test_legacy_umap_endpoint(lib):
    # /molecules/umap still refers to columns from the first (phosphine) database.
    assert lib.get("/molecules/umap").status_code == 200


# --- Identifiers: with and without a pubchem_id column -----------------------------

def test_identifiers(lib):
    r = lib.get("/molecules/identifiers/", smiles=ARYL_4_STORED)
    assert r.status_code == 200
    (body,) = r.json()
    assert body["smiles"] == ARYL_4_CANONICAL
    assert body["InChI"] == "InChI=1S/C6H4N2/c7-5-6-1-3-8-4-2-6/h1-4H"
    assert body["InChIKey"] == ARYL_4_INCHIKEY


def test_identifiers_without_pubchem_column(demo):
    # No pubchem_id column: the backend asks PubChem, and gives None if it cannot
    # reach it. Either way the request succeeds.
    (body,) = demo.get("/molecules/identifiers/", smiles=ARYL_4_STORED).json()
    assert body["PubChemCID"] in (None, 7506)


def test_identifiers_cid_from_database(numeric):
    # The fixture stores a made-up CID for molecule 1, so this value can only
    # have come from the pubchem_id column.
    (body,) = numeric.get("/molecules/identifiers/", smiles=ARYL_1_CANONICAL).json()
    assert body["PubChemCID"] == 999000001


def test_identifiers_empty_pubchem_column(numeric):
    # Column present but empty for this molecule: falls back to PubChem.
    r = numeric.get("/molecules/identifiers/", smiles=ARYL_4_STORED)
    assert r.status_code == 200
    assert r.json()[0]["PubChemCID"] in (None, 7506)


def test_identifiers_invalid_smiles(lib):
    assert lib.get("/molecules/identifiers/", smiles="bad_smiles").status_code == 400


# --- Conformers -----------------------------------------------------------------------

def test_conformer_data(lib):
    r = lib.get("/conformers/data/1")
    assert r.status_code == 200
    body = r.json()
    assert body["molecule_id"] == lib.id(1)
    assert body["method"] == "B3LYP"
    assert "coords" not in body and "elements" not in body


@pytest.mark.parametrize("path", ["/conformers/export/xyz/1", "/conformers/export/1.xyz"])
def test_conformer_xyz(lib, path):
    r = lib.get(path)
    assert r.status_code == 200
    lines = r.text.splitlines()
    assert lines[0] == "15"
    assert len(lines) == 17
    assert lines[2].split()[0] == "N"


@pytest.mark.parametrize("format", ["sdf", "mol"])
def test_conformer_molfile_formats(lib, format):
    r = lib.get(f"/conformers/export/1.{format}")
    assert r.status_code == 200
    assert "V2000" in r.text


@pytest.mark.xfail(strict=True, reason="utils.obabel always writes SDF, whatever format is asked for")
def test_conformer_pdb(lib):
    assert "HETATM" in lib.get("/conformers/export/1.pdb").text


def test_conformer_unknown_format(lib):
    assert lib.get("/conformers/export/1.abc").text == "Not implemented yet"


@known_500
def test_missing_conformer_is_404(lib):
    assert lib.get("/conformers/data/999").status_code == 404


@known_500
def test_other_conformers(lib):
    r = lib.get("/conformers/others_id/1")
    assert r.status_code == 200
    assert r.json() == [1]
