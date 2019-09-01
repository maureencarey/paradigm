#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --time=4-00:00:00
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6000
#SBATCH --account=csblrivanna
#SBATCH --output=/home/mac9jc/paradigm/update_universal_r1b.out

module load anaconda/5.2.0-py3.6
source activate paradigm_env
#python3 ./data_acquistion/update_universal_reaction_set1a.py
python3 ./data_acquistion/update_universal_reaction_set1b.py


