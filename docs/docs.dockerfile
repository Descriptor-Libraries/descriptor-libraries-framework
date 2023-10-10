FROM python:3.11

WORKDIR /docs/
ADD . /docs/

RUN apt-get update \
    && apt-get install -y inotify-tools \
    && pip install -r requirements.txt