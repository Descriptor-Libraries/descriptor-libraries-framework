FROM python

WORKDIR /docs/
ADD . /docs/

RUN apt-get update \
    && apt-get install -y inotify-tools \
    && pip install -r requirements.txt