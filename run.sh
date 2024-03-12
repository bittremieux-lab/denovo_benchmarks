#!/bin/bash
echo "Spectra data dir: $1"
ls $1*.mgf

echo "BUILD ALGORITHM DOCKER"
docker build -f algorithms/casanovo/Dockerfile -t casanovo-docker .

echo "RUN ALGORITHM DOCKER"
docker run --name casanovo_container casanovo_docker $1*.mgf

echo "EXPORT PREDICTIONS"
docker cp casanovo_container:/app/outputs.csv ./

echo "EVALUATE PREDICTIONS"
python3 evaluate.py outputs.csv $1
