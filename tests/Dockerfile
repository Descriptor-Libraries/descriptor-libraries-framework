# Small image that runs the API tests against the backends started by
# docker-compose.test.yml.
FROM python:3.12-slim
RUN pip install --no-cache-dir pytest==8.4.2 requests==2.32.5
WORKDIR /tests
COPY conftest.py test_backend.py ./
CMD ["pytest", "-v"]
