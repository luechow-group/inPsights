# Creating a Docker image for GitLabl CI

Open the directory containing the ```Dockerfile``` and execute the following commands
```
docker build -t registry.git.rwth-aachen.de/luechow-group/amolqcpp .
docker push registry.git.rwth-aachen.de/luechow-group/amolqcpp
```

The container will end up in the
[GitLab Container Registry][https://git.rwth-aachen.de/luechow-group/Amolqcpp/container_registry]
and is accessible in the ```.gitlab-ci.yml``` file with

```
image: registry.git.rwth-aachen.de/luechow-group/amolqcpp:latest
```