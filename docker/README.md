# Metagenome and metatranscriptome sequencing analyses docker

## Dockerfile configuration
The docker image for this workflow is based on the (micromamba docker)[https://github.com/mamba-org/micromamba-docker]. Two micromamba environments exist on this image one containing metawrap 1.2 (the `metawrap_env` environment) and its dependencies, and another (the `metagenome` environment) containing everything else. See the conda enviroment specification files in the asset directory for more detail.

## Execution context switching of the docker image
Micromamba makes use of the environmental variable `ENV_NAME` to determine which environment to activate on entry to the docker container. The default of this variable points to the `metagenome` environment, if you specified the environmental variable `ENV_NAME` as `metawrap_env` you would get the metawrap execution context. 

Alternatively, micromamba can be invoked with context directly, i.e.:
<pre>
micromamba run -n metagenome drep <i>args</i>
micromamba run -n metawrap_env metawrap <i>args</i>
</pre>

## Steps for Building Docker Images
Directly below are instructions for building an image using the provided Dockerfile:

```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker build --no-cache -f example.dockerfile --tag=example:v0.1.0 .

# Testing, take a peek inside
docker run -ti example:v0.1.0 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag example:v0.1.0 skchronicles/example:v0.1.0
docker tag example:v0.1.0 skchronicles/example         # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push skchronicles/example:v0.1.0
docker push skchronicles/example:latest
```

### Other Recommended Steps

Scan your image for known vulnerabilities:

```bash
docker scan example:v0.1.0
```

> **Please Note**: Any references to `skchronicles` should be replaced your username if you would also like to push the image to a non-org account.
