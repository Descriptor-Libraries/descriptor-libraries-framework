# Phosphines WebApp

This web app is made up of several Docker containers which work together through Docker Compose.

## Docker Containers

Images which I have working are marked in <span style="color:green">green</span>, while images which I have not been able to build are marked in <span style="color:red">red.</span> Changes are marked in **bold**.

1.: red_circle: nginx :red_circle: : This is the front end of the app. I cannot build this currently because of memory limitations.

2. <span style="color:green">database</span>: This contains the database connection. **I changed the base image for this image.** It as originally built from a file which was in `backend/database.dockerfile`, with a postgres base image. RDKit was then built and installed. I switched to using a postgres-rdkit base image so that rdkit does not have to be built and installed.

3. <span style="color:green">backend</span>: This is the container which crates the REST API. Uses sqlalchemy, pydantic, and fast api. **I switched to a miniconda base image for this (from fastapi base image).** Allows to more easily install and update Python packages.

4. <span style>cdk-depict</span>: I think something with molecule rendering. I think this is working.

## More about Docker Set-Up

### Volumes
I have set up two volumes to run with these docker containers. These are accessible locally.

1. `kraken-postgres`: Contains database. This database has been updated to contain `umap` information.

2. `pr3`: Contains dumped data (file the database was constructed from)

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
I have added the `umap` information to the database using the [POINT](https://www.postgresql.org/docs/current/datatype-geometric.html) data type. This will allow us to search and filter by distance.

### MolSSI Changes to Code

