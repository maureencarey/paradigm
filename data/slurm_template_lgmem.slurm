#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --time=2-00:00:00
#SBATCH --partition=largemem
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20000
#SBATCH --account=tumi
#SBATCH --output=/home/mac9jc/paradigm/data/slurm_scripts_for_each_species/%x-%j.out

module purge
module load anaconda/5.2.0-py3.6
source activate paradigm_env

