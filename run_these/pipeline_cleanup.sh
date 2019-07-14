#!/bin/bash
# home/mac9jc/paradigm

# CHECK ALL LOGS FOR INFEASIBLE STEPS
cd ~/paradigm/model_generation_logs # paradigm directory
for filename in step*.log; do
echo "$filename"
python3 ../slurm_scripts/check_logs.py $filename
done
