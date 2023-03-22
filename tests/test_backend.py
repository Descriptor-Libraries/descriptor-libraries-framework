"""
Tests for backend
"""
from unicodedata import category
import urllib.parse
import requests
import pytest

# Helper functions to reuse API requests
def get_molecule(molecule_id):

    response = requests.get(f"http://localhost/api/molecules/{molecule_id}")
    return response


def get_molecule_umap(category, show_ml, limit=1000):

    response = requests.get(f"http://localhost/api/molecules/umap?limit={limit}&category={category}&show_ml={show_ml}")
    return response


def get_molecule_substructure_search(substructure, skip=0, limit=100):

    response = requests.get(f"http://localhost/api/molecules/search/?substructure={substructure}&skip={skip}&limit={limit}")
    return response


def get_molecule_neighbors(molecule_id, type, components, skip=0, limit=100):

    response = requests.get(f"http://localhost/api/molecules/{molecule_id}/neighbors/?type={type}&components={components}&skip={skip}&limit={limit}")
    return response


def get_dimensions(type, category, components, skip=0, limit=100):

    response = requests.get(f"http://localhost/api/molecules/dimensions/?type={type}&components={components}&category={category}&skip={skip}&limit={limit}")
    return response


# Tests
# Retrieve tests (status_code == 200)
@pytest.mark.parametrize("molecule_id", [("1"), ("50"), ("12000"), ("234023"), ("331422")])
def test_retrieve_molecule(molecule_id):

    response = get_molecule(molecule_id)
    assert response.status_code == 200


# Note: Wrong categories are returned as an empty JSON.
@pytest.mark.parametrize("category, show_ml", [("", False), ("", True), ("pc3", False), ("pcn", False), ("pon", True), ("wrong_category", True)])
def test_retrieve_molecule_umap(category, show_ml):

    response = get_molecule_umap(category, show_ml)
    assert response.status_code == 200


@pytest.mark.parametrize("substructure", [("C=O"), ("C=CC"), ("C=N"), ("N=N"), ("C-P"), ("C=C-P")])
def test_retrieve_substructure_search(substructure):

    response = get_molecule_substructure_search(substructure)
    assert response.status_code == 200


@pytest.mark.parametrize("molecule_id, type, components", [(1, "pca", "1, 2, 3"), (1, "umap", "1, 2"), (1, "pca", "1, 3, 2"), (1, "pca", "1, 4"), (1, "pca", "1"), (1, "umap", "1"), (35, "umap", "1, 2"), (300000, "umap", "2"), (35, "PCA", "1, 2"), (35, "UMAP", "1, 2")])
def test_retrieve_molecule_neighbors(molecule_id, type, components):

    response = get_molecule_neighbors(molecule_id, type, components)
    assert response.status_code == 200


@pytest.mark.parametrize("type, category, components", [("pca", "", "1, 2, 3"), ("pca", "", "1"), ("pca", "", "4"), ("pca", "", "1, 3, 2"), ("pca", "", "1, 4"), ("umap", "", "1, 2"), ("umap", "", "1"), ("umap", "", "2"), ("umap", "pc3", "1, 2"), ("umap", "pon", "1, 2"), ("umap", "pcn", "1, 2")])
def test_retrieve_dimensions(type, category, components):

    response = get_dimensions(type, category, components)
    assert response.status_code == 200


# TODO:
# Failing tests (status_code == 400)
# TODO: Set up catch on backend if the molecule id is out of range.
@pytest.mark.parametrize("molecule_id", [("0"), ("-1"), ("-50"), ("-12000"), ("-234023"), ("-331422")])
def test_retrieve_molecule_under_range(molecule_id):

    response = get_molecule(molecule_id)
    assert response.status_code == 500


# For numbers greater than we have it should be a 404.
@pytest.mark.parametrize("molecule_id", [("331423"), ("350000"), ("400000")])
def test_retrieve_molecule_over_range(molecule_id):

    response = get_molecule(molecule_id)
    assert response.status_code == 500


@pytest.mark.parametrize("substructure", [("C=XX"), ("P=XX"), ("P=F"), ("C==C"), ("-C-"), ("C=F=C")])
def test_retrieve_substructure_search_invalid_smiles(substructure):

    response = get_molecule_substructure_search(substructure)
    assert response.status_code == 400


@pytest.mark.parametrize("molecule_id, type, components", [(1, "wrong_type", "1, 2, 3"), (1, "pcca", "1, 2, 3"), (1, "umaaaap", "1, 2, 3"), (1, "pca", "1, 2, 3, 4, 5"), (1, "pca", "1, 15"), (1, "pca", "45"), (1, "umap", "3")])
def test_retrieve_molecule_neighbors_invalid_parameters(molecule_id, type, components):

    response = get_molecule_neighbors(molecule_id, type, components)
    assert response.status_code == 400


@pytest.mark.parametrize("type, category, components", [("wrong_type", "", "1, 2, 3"), ("pca", "wrong_category", "1, 2, 3"), ("pca", "", "1, 2, 3, 4, 5"), ("pca", "", "1, 15"), ("pca", "", "45"), ("umap", "", "3")])
def test_retrieve_dimensions_invalid_parameters(type, category, components):

    response = get_dimensions(type, category, components)
    assert response.status_code == 400


# Output tests (correct number of items returned)
# Output tests (correct components returned)
# Output tests (correct values returned)

