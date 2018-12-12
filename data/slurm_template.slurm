#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --time=2-00:00:00
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6000
#SBATCH --account=csblrivanna

module load anaconda3
module load gurobi/8.0.1
source activate paradigm_env

