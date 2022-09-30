FROM continuumio/miniconda3 as react-build

RUN yes | conda install -c conda-forge nodejs">17"
RUN yes | conda install -c conda-forge yarn

WORKDIR /app/
ADD ./frontend /app/
EXPOSE 3000
RUN yarn add react-scripts
RUN yarn
