## Steps for Building Docker Images

Directly below are instructions for building an image using the provided Dockerfile:

```bash
# See listing of images on computer
docker image ls

# Build from Dockerfile
docker buildx build --platform linux/amd64 --no-cache -f Dockerfile --tag=metaphlan-utils:v0.1.0 .

# Testing, take a peek inside
docker run --platform linux/amd64 -ti metaphlan-utils:v0.1.0 /bin/bash

# Updating Tag  before pushing to DockerHub
docker tag metaphlan-utils:v0.1.0 skchronicles/metaphlan-utils:v0.1.0
docker tag metaphlan-utils:v0.1.0 skchronicles/metaphlan-utils         # latest

# Check out new tag(s)
docker image ls

# Push new tagged image to DockerHub
docker push --platform linux/amd64 skchronicles/metaphlan-utils:v0.1.0
docker push --platform linux/amd64 skchronicles/metaphlan-utils:latest
```

### Other Recommended Steps

Scan your image for known vulnerabilities:

```bash
docker scan metaphlan-utils:v0.1.0
```

> **Please Note**: Any references to `skchronicles` should be replaced your username if you would also like to push the image to a registry like DockerHub.