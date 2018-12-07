#!/bin/bash

module load gurobi-7.0.2
module load anaconda3
source activate gurobienv
cd ./paradigm/data

mkdir ./model_generation_logs
mv ./*.log ./model_generation_logs
mkdir ./models
mv ./*.json ./models

# check logs for infeasible steps
cd ./model_generation_logs
for filename in *2018.log; do
    echo "$filename"
    python ../check_logs.py $filename
done
