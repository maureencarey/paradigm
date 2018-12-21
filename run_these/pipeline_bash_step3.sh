#!/bin/bash
# path = "/home/mac9jc/paradigm"

cd /home/mac9jc/paradigm/data

mkdir /home/mac9jc/paradigm/model_generation_logs
mv ./*.log /home/mac9jc/paradigm/model_generation_logs
mkdir ./models
mv ./*.json ./models

# check logs for infeasible steps
cd /home/mac9jc/paradigm/model_generation_logs
for filename in *2018.log; do
    echo "$filename"
    python ../slurm_scripts/check_logs.py $filename
done

#git add .
#git commit -m "ran full pipeline on Dec 21th, 2018"
#git push origin master

