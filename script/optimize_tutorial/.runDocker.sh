#!/bin/bash

docker build -t qdpy_docker . --no-cache
# docker build -t qdpy_docker . 
currentdir=$(pwd)
basename=`basename $currentdir`
docker run  --rm -it -v "${currentdir}:/home/user/${basename}" qdpy_docker 
# 