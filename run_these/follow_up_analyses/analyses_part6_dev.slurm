#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mac9jc@virginia.edu
#SBATCH --partition=dev
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6000
#SBATCH --account=csblrivanna
#SBATCH --output=%x-%j.out

module purge
module load anaconda/5.2.0-py3.6
source activate paradigm_env
python3 ./model_evaluation/biomass_precursors.py
