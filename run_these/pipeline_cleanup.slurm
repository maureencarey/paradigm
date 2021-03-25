#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mac9jc@virginia.edu
#SBATCH --time=2-00:00:00
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=9000
#SBATCH --account=csblrivanna
#SBATCH --output=/home/mac9jc/paradigm/generating_check_log_and_outputs.out

module purge
module load anaconda/5.2.0-py3.6
source activate paradigm_env

cd ~/paradigm/model_generation_logs # paradigm directory
for filename in step*.log; do
#    echo "$filename"
    python3 ../slurm_scripts/check_logs.py $filename
done

cd ~/paradigm/data/slurm_outputs
for filename in *.out; do
#    echo "$filename"
    python3 ~/paradigm/slurm_scripts/check_outputs.py $filename
done

echo 'done'
