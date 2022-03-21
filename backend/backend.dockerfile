FROM tiangolo/uvicorn-gunicorn-fastapi:python3.7

WORKDIR /app/

COPY ./requirements.txt requirements.txt
RUN pip install -r requirements.txt && rm -fr requirements.txt

COPY ./app /app
ENV PYTHONPATH=/app
