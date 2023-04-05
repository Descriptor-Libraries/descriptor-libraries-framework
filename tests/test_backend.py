"""
Tests for backend
"""
import requests
import pytest
import json


# Helper functions to reuse API requests
def get_molecule(molecule_id):
    response = requests.get(f"http://localhost/api/molecules/{molecule_id}")
    return response


def get_molecule_umap(category, show_ml, limit=1000):
    response = requests.get(
        f"http://localhost/api/molecules/umap?limit={limit}&category={category}&show_ml={show_ml}"
    )
    return response


def get_molecule_substructure_search(substructure, skip=0, limit=100):
    response = requests.get(
        f"http://localhost/api/molecules/search/?substructure={substructure}&skip={skip}&limit={limit}"
    )
    return response


def get_molecule_neighbors(molecule_id, type, components, skip=0, limit=100):
    response = requests.get(
        f"http://localhost/api/molecules/{molecule_id}/neighbors/?type={type}&components={components}&skip={skip}&limit={limit}"
    )
    return response


def get_dimensions(type, category, components, skip=0, limit=100):
    response = requests.get(
        f"http://localhost/api/molecules/dimensions/?type={type}&components={components}&category={category}&skip={skip}&limit={limit}"
    )
    return response


# Tests
# Retrieve tests (status_code == 200)
@pytest.mark.parametrize("molecule_id", [("1"), ("129"), ("1290"), ("57"), ("565")])
def test_retrieve_molecule(molecule_id):
    response = get_molecule(molecule_id)
    assert response.status_code == 200


# Note: Wrong categories are returned as an empty JSON.
@pytest.mark.parametrize(
    "category, show_ml",
    [
        ("", False),
        ("", True),
        ("pc3", False),
        ("pcn", False),
        ("pon", True),
        ("wrong_category", True),
    ],
)
def test_retrieve_molecule_umap(category, show_ml):
    response = get_molecule_umap(category, show_ml)
    assert response.status_code == 200


@pytest.mark.parametrize(
    "substructure", [("C=O"), ("C=CC"), ("C=N"), ("N=N"), ("C-P"), ("C=C-P")]
)
def test_retrieve_substructure_search(substructure):
    response = get_molecule_substructure_search(substructure)
    assert response.status_code == 200


@pytest.mark.parametrize(
    "molecule_id, type, components",
    [
        (1, "pca", "1, 2, 3"),
        (1, "umap", "1, 2"),
        (1, "pca", "1, 3, 2"),
        (1, "pca", "1, 4"),
        (1, "pca", "1"),
        (1, "umap", "1"),
        (57, "umap", "1, 2"),
        (565, "umap", "2"),
        (57, "PCA", "1, 2"),
        (57, "UMAP", "1, 2"),
    ],
)
def test_retrieve_molecule_neighbors(molecule_id, type, components):
    response = get_molecule_neighbors(molecule_id, type, components)
    assert response.status_code == 200


@pytest.mark.parametrize(
    "type, category, components",
    [
        ("pca", "", "1, 2, 3"),
        ("pca", "", "1"),
        ("pca", "", "4"),
        ("pca", "", "1, 3, 2"),
        ("pca", "", "1, 4"),
        ("umap", "", "1, 2"),
        ("umap", "", "1"),
        ("umap", "", "2"),
        ("umap", "pc3", "1, 2"),
        ("umap", "pcn", "1, 2"),
    ],
)
def test_retrieve_dimensions(type, category, components):
    response = get_dimensions(type, category, components)
    assert response.status_code == 200


# Failing tests (status_code == 400)
@pytest.mark.parametrize(
    "molecule_id", [("0"), ("-1"), ("-50"), ("-12000"), ("-234023"), ("-331422")]
)
def test_retrieve_molecule_under_range(molecule_id):
    response = get_molecule(molecule_id)
    assert response.status_code == 500


# For numbers greater than we have it should be a 404.
@pytest.mark.parametrize("molecule_id", [("331423"), ("350000"), ("400000")])
def test_retrieve_molecule_over_range(molecule_id):
    response = get_molecule(molecule_id)
    assert response.status_code == 404


@pytest.mark.parametrize(
    "substructure", [("C=XX"), ("P=XX"), ("P=F"), ("C==C"), ("-C-"), ("C=F=C")]
)
def test_retrieve_substructure_search_invalid_smiles(substructure):
    response = get_molecule_substructure_search(substructure)
    assert response.status_code == 400


@pytest.mark.parametrize(
    "molecule_id, type, components",
    [
        (1, "wrong_type", "1, 2, 3"),
        (1, "pcca", "1, 2, 3"),
        (1, "umaaaap", "1, 2, 3"),
        (1, "pca", "1, 2, 3, 4, 5"),
        (1, "pca", "1, 15"),
        (1, "pca", "45"),
        (1, "umap", "3"),
    ],
)
def test_retrieve_molecule_neighbors_invalid_parameters(molecule_id, type, components):
    response = get_molecule_neighbors(molecule_id, type, components)
    assert response.status_code == 400


@pytest.mark.parametrize(
    "type, category, components",
    [
        ("wrong_type", "", "1, 2, 3"),
        ("pca", "wrong_category", "1, 2, 3"),
        ("pca", "", "1, 2, 3, 4, 5"),
        ("pca", "", "1, 15"),
        ("pca", "", "45"),
        ("umap", "", "3"),
        ("umap", "pon", "1, 2"),
    ],
)
def test_retrieve_dimensions_invalid_parameters(type, category, components):
    response = get_dimensions(type, category, components)
    assert response.status_code == 400


# Output tests (correct values returned)
@pytest.mark.parametrize("molecule_id", [("1"), ("129"), ("1290"), ("57"), ("565")])
def test_retrieve_molecule_output_results(molecule_id):
    response = get_molecule(molecule_id)

    # Exporting JSONs once to have ground truth, only generate once!
    # with open(f'data/get_molecule_{molecule_id}.json', 'w', encoding='utf-8') as f:
    #     json.dump(response.json(), f, ensure_ascii=False, indent=4)

    # Read ground truth data and then compare.
    with open(f"data/get_molecule_{molecule_id}.json") as json_file:
        data = json.load(json_file)

    # Check to see if the data returned is identical
    assert response.json() == data


@pytest.mark.parametrize(
    "test_case, category, show_ml",
    [
        (1, "", False),
        (2, "", True),
        (3, "pc3", False),
        (4, "pcn", False),
        (5, "pon", True),
        (6, "wrong_category", True),
    ],
)
def test_retrieve_molecule_umap_results(test_case, category, show_ml):
    response = get_molecule_umap(category, show_ml)

    # Exporting JSONs once to have ground truth, only generate once!
    # with open(f'data/get_molecule_umap_{test_case}.json', 'w', encoding='utf-8') as f:
    #     json.dump(response.json(), f, ensure_ascii=False, indent=4)

    # Read ground truth data and then compare.
    with open(f"data/get_molecule_umap_{test_case}.json") as json_file:
        data = json.load(json_file)

    # Check to see if the data returned is identical
    assert response.json() == data


@pytest.mark.parametrize(
    "test_case, substructure",
    [(1, "C=O"), (2, "C=CC"), (3, "C=N"), (4, "N=N"), (5, "C-P"), (6, "C=C-P")],
)
def test_retrieve_substructure_search_results(test_case, substructure):
    response = get_molecule_substructure_search(substructure)

    # Exporting JSONs once to have ground truth, only generate once!
    # with open(f'data/get_molecule_substructure_search_{test_case}.json', 'w', encoding='utf-8') as f:
    #     json.dump(response.json(), f, ensure_ascii=False, indent=4)

    # Read ground truth data and then compare.
    with open(f"data/get_molecule_substructure_search_{test_case}.json") as json_file:
        data = json.load(json_file)

    # Check to see if the data returned is identical
    assert response.json() == data


@pytest.mark.parametrize(
    "test_case, molecule_id, type, components",
    [
        (1, 1, "pca", "1, 2, 3"),
        (2, 1, "umap", "1, 2"),
        (3, 1, "pca", "1, 3, 2"),
        (4, 1, "pca", "1, 4"),
        (5, 1, "pca", "1"),
        (6, 57, "umap", "1"),
        (7, 57, "umap", "1, 2"),
        (8, 565, "umap", "2"),
        (9, 57, "PCA", "1, 2"),
        (10, 57, "UMAP", "1, 2"),
    ],
)
def test_retrieve_molecule_neighbors_results(test_case, molecule_id, type, components):
    response = get_molecule_neighbors(molecule_id, type, components)

    # Exporting JSONs once to have ground truth, only generate once!
    # with open(f'data/get_molecule_neighbors_{test_case}.json', 'w', encoding='utf-8') as f:
    #     json.dump(response.json(), f, ensure_ascii=False, indent=4)

    # Read ground truth data and then compare.
    with open(f"data/get_molecule_neighbors_{test_case}.json") as json_file:
        data = json.load(json_file)

    # Check to see if the data returned is identical
    assert response.json() == data


@pytest.mark.parametrize(
    "test_case, type, category, components",
    [
        (1, "pca", "", "1, 2, 3"),
        (2, "pca", "", "1"),
        (3, "pca", "", "4"),
        (4, "pca", "", "1, 3, 2"),
        (5, "pca", "", "1, 4"),
        (6, "umap", "", "1, 2"),
        (7, "umap", "", "1"),
        (8, "umap", "", "2"),
        (9, "umap", "pc3", "1, 2"),
        (10, "pca", "pcn", "3,1"),
        (11, "umap", "pcn", "1, 2"),
    ],
)
def test_retrieve_dimensions_results(test_case, type, category, components):
    response = get_dimensions(type, category, components)

    # Exporting JSONs once to have ground truth, only generate once!
    # with open(f'data/get_dimensions_{test_case}.json', 'w', encoding='utf-8') as f:
    #     json.dump(response.json(), f, ensure_ascii=False, indent=4)

    # Read ground truth data and then compare.
    with open(f"data/get_dimensions_{test_case}.json") as json_file:
        data = json.load(json_file)

    # Check to see if the data returned is identical
    assert response.json() == data
