#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mac9jc@virginia.edu
#SBATCH --time=4-00:00:00
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6000
#SBATCH --account=csblrivanna
#SBATCH --output=/home/mac9jc/paradigm/step2_update_universal_r1.out

module purge
module load anaconda/5.2.0-py3.6
conda activate annotation_env
#source activate paradigm_env
python3 ./data_acquistion/update_universal_reaction_set1a.py
python3 ./data_acquistion/update_universal_reaction_set1b.py


