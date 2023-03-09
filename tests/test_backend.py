"""
Tests for backend
"""

import requests

def test_retrieve_molecule():

    response = requests.get("http://localhost/api/molecules/1")
    assert response.status_code == 200