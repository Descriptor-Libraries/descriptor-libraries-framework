"""
Tests for backend
"""
import urllib.parse
import requests

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
def test_retrieve_molecule():

    response = get_molecule(1)
    assert response.status_code == 200


def test_retrieve_molecule_umap():

    response = get_molecule_umap("", False)
    assert response.status_code == 200


def test_retrieve_substructure_search():

    response = get_molecule_substructure_search("C=O")
    assert response.status_code == 200


def test_retrieve_molecule_neighbors():

    response = get_molecule_neighbors(1, "pca", "1, 2, 3")
    assert response.status_code == 200


def test_retrieve_dimensions():

    response = get_dimensions("pca", "", "1, 2, 3")
    assert response.status_code == 200


# TODO:
# Failing tests (status_code == 400)
# Output tests (correct number of items returned)
# Output tests (correct components returned)
# Output tests (correct values returned)

