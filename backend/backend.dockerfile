FROM continuumio/miniconda3

RUN yes | conda install -c conda-forge rdkit fastapi alembic psycopg2-binary sqlalchemy tenacity uvicorn curl

WORKDIR /app/
ADD ./app /app/
RUN pip install -e .

