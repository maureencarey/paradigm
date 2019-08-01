#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --partition=dev
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=9000
#SBATCH --account=csblrivanna

module load anaconda/5.2.0-py3.6
source activate paradigm_env

python3 data_acquistion/update_universal_reaction_set.py universal_model_updated.json
