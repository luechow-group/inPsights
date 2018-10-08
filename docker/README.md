# Creating a Docker image for GitLabl CI

## Prerequisites
Install Docker Community Edition (CE)on your local machine via
[command line](https://docs.docker.com/install/linux/docker-ce/ubuntu/#install-docker-ce-1)
or
[download it as a binary](https://store.docker.com/search?type=edition&offering=community).

## Login to the GitLab Docker registry
Start the docker daemon and login to the GitLab Docker registry with your GitLab credentials
(GitLab username and password):
```
docker login registry.git.rwth-aachen.de
```

## Build and push to the GitLab Docker registry
Open the directory containing the ```Dockerfile``` and execute the following commands
```
docker build -t registry.git.rwth-aachen.de/luechow-group/amolqcpp .
docker push registry.git.rwth-aachen.de/luechow-group/amolqcpp
```

The container will show up in the
[GitLab Container Registry](https://git.rwth-aachen.de/luechow-group/Amolqcpp/container_registry)
and is accessible in the ```.gitlab-ci.yml``` file with
```
image: registry.git.rwth-aachen.de/luechow-group/:latest
```
