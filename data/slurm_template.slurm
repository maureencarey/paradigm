#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --time=2-00:00:00
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=12000
#SBATCH --account=csblrivanna
#SBATCH --output=%x-%j.out

module purge
module load anaconda/5.2.0-py3.6
conda activate annotation_env
#source activate paradigm_env

