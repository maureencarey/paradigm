#!/bin/bash
# home/mac9jc/paradigm

# cd home/mac9jc/paradigm
#mkdir home/mac9jc/paradigm/model_generation_logs
#mv ./data/*.log ./model_generation_logs

# mkdir ./models
#mv ./data/*.json ./models
#mv ./data/*.xml ./models

# CHECK ALL LOGS FOR INFEASIBLE STEPS
cd ~/paradigm/model_generation_logs # paradigm directory
for filename in step*.log; do
echo "$filename"
python3 ../slurm_scripts/check_logs.py $filename
done
