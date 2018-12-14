#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --time=2-00:00:00
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6000
#SBATCH --account=csblrivanna

module load R/3.4.3
module load diamond
bash ./slurm_scripts/pipeline_slurm_step1a.sh

module load anaconda3
source activate paradigm_env
module load gurobi/8.0.1
bash ./slurm_scripts/pipeline_slurm_step1a.sh
