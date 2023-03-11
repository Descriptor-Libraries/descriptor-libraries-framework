Architecture Overview
=====================

kraken is an app consisting of many different services stitched together with `Docker Compose <https://docs.docker.com/compose/>`_.

The site consists of the following services:

#. A pre-built container for postgresql database with an RDKit plugin.
#. A Python "backend" that serves a REST API using FastAPI.
#. A React front end for the graphical portion of the website
#. A container running Sphinx for documentation.
#. A pre-built container, `cdk-depict <https://github.com/cdk/depict>`_, which adds an API for retrieving depicitions of molecules based on SMILES strings. 
#. A reverse proxy to control routing, and make all the services appear as part of the same site.
