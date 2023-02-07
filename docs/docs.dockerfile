FROM python

WORKDIR /docs/
ADD . /docs/

RUN apt update
RUN yes | apt install inotify-tools
RUN yes | pip install -r requirements.txt