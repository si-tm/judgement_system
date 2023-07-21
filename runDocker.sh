#!/bin/bash

# docker build -t qdpy_docker . --no-cache
docker build -t qdpy_docker . 
currentdir=$(pwd)
basename=`basename $currentdir`
# docker run  --rm -it -v "${currentdir}:/home/user/${basename}" qdpy_docker 
# docker run  --rm -v "${currentdir}:/home/user/${basename}" qdpy_docker 
docker run  --rm -it -v "${currentdir}:/home/user/${basename}" qdpy_docker 
# cd /home/user/judgement_system/script/optimize_tutorial
