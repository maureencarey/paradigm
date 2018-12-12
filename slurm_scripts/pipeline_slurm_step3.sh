#!/bin/bash

module load gurobi/8.0.1
module load anaconda3
source activate paradigm_env

cd ./paradigm/data

mkdir /scratch/mac9jc/paradigm/model_generation_logs
mv ./*.log /scratch/mac9jc/paradigm/model_generation_logs
mkdir ./models
mv ./*.json ./models

# check logs for infeasible steps
cd /scratch/mac9jc/paradigm/model_generation_logs
for filename in *2018.log; do
    echo "$filename"
    python ../slurm_scripts/check_logs.py $filename
done

git add .
git commit -m "ran full pipeline on Dec 12th, 2018"
git push origin master
