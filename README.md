# Phosphines WebApp

This web app is made up of several Docker containers which work together through Docker Compose.

## Docker Containers

Images which I have working are marked with a green circle :green_circle:, while images which I have not been able to build are marked with a red circle :red_circle". Changes are marked in **bold**.

1. :red_circle: nginx :red_circle: : This is the front end of the app. I cannot build this currently because of memory limitations.

2. :green_circle: database :green_circle: : This contains the database connection. **I changed the base image for this image.** It as originally built from a file which was in `backend/database.dockerfile`, with a postgres base image. RDKit was then built and installed. I switched to using a postgres-rdkit base image so that rdkit does not have to be built and installed.

3. :green_circle: backend :green_circle: This is the container which crates the REST API. Uses sqlalchemy, pydantic, and fast api. **I switched to a miniconda base image for this (from fastapi base image).** Allows to more easily install and update Python packages.

4. <span style>cdk-depict</span>: I think something with molecule rendering. I think this is working.

## More about Docker Set-Up

### Volumes
I have set up three volumes to run with these docker containers. The first two are accessible locally on my set-up, but are not included in the remote repository. The third volume mounts the Python code so it can be developed while the containers are run. This will likely be removed in production.

1. `kraken-postgres`: Contains database. This database has been updated to contain `umap` information.

2. `pr3`: Contains dumped data (file the database was constructed from using `pg_restore`)

3. backend volume for Python code to allow development while container is running.

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

To stop the containers, use:

```
docker-compose down
```

## Working with Postgres Database
Once the database docker container is running, you can connect to the database from outside of the container using `psql`, if installed.

```
psql -U postgres -h localhost
```
Connecting to the database from outside of the Docker container will require entering the admin password.

Alternatively, the database docker container can be entered, and a password will not be required to work with the database. Use `docker ps` to get the container ID of the database container, then you can do

```
docker exec -it CONTAINER_ID /bin/bash
```

To start an interactive bash terminal in the running container.

### Database
1. This database was created using a database dump file and supplemented with PCA and UMAP data provided in CSV format. This docker compose file assumes that a database with relevant data already exists and is mounted on a volume called `kraken-postgres`.

#### Importing from CSV
The approach I've taken is to use pandas to create a csv which has `smiles` and the data I want to import into the database. 

Then, I log into the database and create a new temporary table. For example, for the `umap` data,

```
CREATE TABLE temporary
(smiles VARCHAR PRIMARY KEY,
umap1 NUMERIC,
umap2 NUMERIC,
umap POINT
);
```

Then, to read in the data from a CSV

```
\copy temporary(smiles, umap1, umap2)
FROM '\PATH\TO\CSV'
DELIMITER ','
CSV HEADER;
```

To create the point type

```
UPDATE new_data
SET umap = point(umap1, umap2);
```

Add column to molecule table

```
ALTER TABLE molecule
ADD COLUMN umap point;
```

Set column equal to `umap` from temporary table.

```
UPDATE molecule
SET umap = nd.map
FROM temporary nd
WHERE nd.smiles = molecule.smiles;
```

### MolSSI Changes to the Database
I have added the `umap` information to the database using the [POINT](https://www.postgresql.org/docs/current/datatype-geometric.html) data type. The point data type is for storing two dimensional data. This will allow us to search and filter by umap distance.

### MolSSI Changes to Code

