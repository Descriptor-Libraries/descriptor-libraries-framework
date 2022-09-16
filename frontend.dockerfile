FROM continuumio/miniconda3

RUN yes | conda install -c conda-forge nodejs">17"
RUN yes | conda install -c conda-forge yarn

RUN conda list nodejs

WORKDIR /app/
ADD ./frontend /app/
EXPOSE 3000
RUN yarn
