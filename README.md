# kraken web app

This repository provides a web application for displaying molecule library information. Our use case is "kraken - Kolossal viRtual dAtabase formoleKular dEscriptors of orgaNophosphorus ligands." You can see the application at [kraken.molssi.org](https://kraken.molssi.org/,)
The app is made up of several Docker containers which work together through Docker Compose. 
The docker containers connect to a database which contains molecule + conformer information.

Running the application locally
================================

Running the kraken app locally requires installation of Docker and Docker Compose.
Follow the instructions from Docker to install [Docker Engine](https://docs.docker.com/engine/install/) 
and install [Docker Compose](https://docs.docker.com/compose/install/>) for your operating system.
Next, clone the [project repository](https://github.com/janash/kraken-starting) to your local machine.

Running the testing app - (self contained)
------------------------------------------
Once you have Docker, Docker Compose, and have cloned the repository, navigate to the top level of the project.
You can run the app using a minimal data set by running the "testing app". 

```
docker-compose -f docker-compose.test.yml up
```
To run the app in the background, you can add the argument `--detach` to the command above. 

The first time you run this, all of the images will be built, and will take a few minutes. After running the command, you should be able to navigate to ``http://localhost`` in your browser and see
the app running. 

To stop the containers, use:

```
docker-compose down
```

## Working with Postgres Database
If you woud like to install a full working development version of the app, you will have to download the full dataset (coming soon), create volumes for the containers, and laod the dataset. 
Running a local version of the full app uses the `docker-compose.local.yml` file

1. Create volumes. 

   ```
   docker volume create kraken-postgres
   docker volume create pr3
   ```
2. Download database dump - coming soon.
3. Start the database container
```
docker-compose -f docker-compose.local.yml run local_database
```
4. Next, you will load the database dump file into the `kraken-4gb` volume. For this, I put the dump file in the pr3 volume, and you can enter the container and [load the dump file](https://www.postgresql.org/docs/current/backup-dump.html
). If you have the postgres sql client installed, you can also do this from the host side (without entering the docker container).

To the database from outside of the container using `psql`, if installed.

```
psql -U postgres -h localhost
```
Connecting to the database from outside of the Docker container will require entering the admin password.

Alternatively, the database docker container can be entered, and a password will not be required to work with the database. Use `docker ps` to get the container ID of the database container, then you can do

```
docker exec -it CONTAINER_ID /bin/bash
```

to start an interactive bash terminal in the running container. From there, you can access the dump file in the pr3 volume (mapped to `home/db`) to load the dump file.

Future Plans
------------
We plan to continue cleaning up this application and implementing new features. Our intention is to eventually create this as a template project "cookiecutter" so that people can easily create applications that display molecular information.

Acknowledgements
----------------
This project is a collaboration between [The Molecular Sciences Software Institute](https://molssi.org/) and the [Center for Computer Assisted Synthesis](https://ccas.nd.edu/).

MolSSI is funded by the National Science Foundation OAC-1547580 and CHE-2136142.
C-CAS is funded by the National Science Foundation CHE–2202693.

