Setting up your development environment
========================================

Running the kraken app locally requires installation of Docker and Docker Compose.
Follow the instructions from Docker to install `Docker Engine <https://docs.docker.com/engine/install/>`_ 
and `install Docker Compose <https://docs.docker.com/compose/install/>`_ for your operating system.
Next, clone the `project repository <https://github.com/janash/kraken-starting>`_ to your local machine.
The project repository is currently private, but should be open and public soon.

Running the testing app
-----------------------
Once you have Docker, Docker Compose, and have cloned the repository, navigate to the top level of the project.
You can run the app using a minimal data set by running the "testing app". 
    
.. code-block:: shell

    docker-compose -f docker-compose.test.yml

To run the app in the background, you can add the argument `--detach` to the command above. 
After running the command, you should be able to navigate to ``http://localhost`` in your browser and see
the app running. The testing app runs the same code as the current "production" ``docker-compose.yml``, but with
a minimal data set. 
You can use this as your development environment if you wish to.
The testing app has a small database with a non-persistent database, 
meaning that if you modify the data in the database, the changes will not be present the next time you start the app.

Running with the full data set
------------------------------
To run with the full data set, you will have to download the full postgres database, 
mount two Docker volumes on your machine, and build the database.

More directions coming soon.