"""
Settings for the backend API tests.

The tests talk to running backends over HTTP. Where they are is a setting:

    API_URL          library with string IDs and no pubchem_id column (the demo data)
    API_URL_NUMERIC  library with integer IDs, a pubchem_id column and ML data

`docker-compose.test.yml` starts both and sets both variables. To run against
the demo stack from the README instead, only the first one is needed:

    API_URL=http://localhost/api/demo pytest tests

Tests that need the second library are skipped when API_URL_NUMERIC is not set.
"""
import os

import pytest
import requests


class Api:
    """A tiny client: `api.get("/molecules/first_id")`."""

    def __init__(self, base_url):
        self.base_url = base_url.rstrip("/")

    def get(self, path, **params):
        return requests.get(self.base_url + path, params=params or None, timeout=30)


@pytest.fixture(scope="session")
def demo():
    return Api(os.environ.get("API_URL", "http://localhost/api/demo"))


@pytest.fixture(scope="session")
def numeric():
    url = os.environ.get("API_URL_NUMERIC")
    if not url:
        pytest.skip("API_URL_NUMERIC is not set")
    return Api(url)
