# Running the project in a docker container

To build and run the Docker container, follow these steps:

- Open a terminal and navigate to the directory where the docker-compose and Dockerfile are saved.
- Build the Docker image by running the following command: 
```
$ docker-compose build
```

Access the container's Bash shell by running the following command

```
$ docker-compose up -d
$ docker-compose exec nextflow-genomic bash
```

Initialise conda

```bash
$ conda init bash
$ source ~/.bashrc
```

Activate the environment
```bash
$ conda activate snakemake-tutorial
```