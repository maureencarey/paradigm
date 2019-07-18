#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --time=2-00:00:00
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=9000
#SBATCH --account=csblrivanna
#SBATCH --output=%x-%j.out

module load anaconda/5.2.0-py3.6
source activate paradigm_env

cd ~/paradigm/model_generation_logs # paradigm directory
for filename in step*.log; do
    echo "$filename"
    python3 ../slurm_scripts/check_logs.py $filename
done

cd ~/paradigm/data/slurm_outputs
for filename in *.out; do
    echo "$filename"
    python3 ../../slurm_scripts/check_outputs.py $filename
done

