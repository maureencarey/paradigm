#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --time=2-00:00:00
#SBATCH --output=/home/mac9jc/paradigm/pipeline_step1.out
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6000
#SBATCH --account=csblrivanna
#SBATCH --output=%x-%j.out

module load gcc/7.1.0
module load R/3.3.0
module load diamond
bash ./slurm_scripts/pipeline_slurm_step1a.sh

module load anaconda/5.2.0-py3.6
source activate paradigm_env
bash ./slurm_scripts/pipeline_slurm_step1b.sh
