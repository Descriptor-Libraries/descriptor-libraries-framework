FROM informaticsmatters/rdkit-python3-debian

WORKDIR /app/

COPY ./requirements.txt requirements.txt
RUN pip install -r requirements.txt

COPY ./app /app
ENV PYTHONPATH=/app
