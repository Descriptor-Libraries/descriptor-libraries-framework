# Phosphines WebApp

This web app is made up of several Docker containers which work together through Docker Compose.

## Docker Containers

Images which I have working are marked in <span style="color:green">green</span>, while images which I have not been able to build are marked in <span style="color:red">red.</span>

1. <span style="color:red">nginx</span>: This is the front end of the app. I cannot build this currently because of memory limitations.

2. <span style="color:green">database</span>: This contains the database connection. I chaned the base image for this image. It as originally built from a file which was in `backend/database.dockerfile`, with a postgres base image. RDKit was then built and installed. I switched to using a postgres-rdkit base image so that rdkit does not have to be built and installed.

3. <span style="color:green">backend</span>: This is the container which crates the REST API. Uses sqlalchemy, pydantic, and fast api. I switched to a miniconda base image for this (from fastapi base image). Allows to more easily install and update Python packages.

4. <span style>cdk-depict</span>: I think something with molecule rendering. I think this is working.

## More about Docker Set-Up

### Volumes
1. `kraken-postgres`: Contains database.

2. `pr3`: Contains dumped data.

### Database
1. This database was created using a database dump file and supplemented with PCA and UMAP data provided in CSV format. This docker compose file assumes that a database with relevant data already exists and is mounted on a volume called `kraken-postgres`.

### Running containers
1. All containers can be built using

```
docker-compose build
```

The argument `--no-cache` can be added to disable the cache during building.

To run and connect all of the containers, use:

```
docker-compose up
```


## Working with Postgres Database
Once the database docker container is running, you can connect to the database from outside of the container using `psql`, if installed.

```
psql -U postgres -h localhost
```
Connecting to the database from outside of the Docker container will require entering the admin password.

Alternatively, the database docker container can be entered, and a password will not be required to work with the database. Use `docker ps` to get the container ID of the database container, then you can do

```
docker exect -it CONTAINER_ID /bin/bash
```

To start an interactive bash terminal in the running container.