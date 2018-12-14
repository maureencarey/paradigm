#!/bin/bash
# home/mac9jc/paradigm

mkdir ./model_generation_logs
mv ./data/*.log ./model_generation_logs
mkdir ./models
mv ./data/*.json ./models

# CHECK ALL LOGS FOR INFEASIBLE STEPS
cd ./model_generation_logs # paradigm directory
for filename in *2018.log; do
echo "$filename"
python ../check_logs.py $filename
done
