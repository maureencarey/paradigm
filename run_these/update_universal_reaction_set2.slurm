#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mac9jc@virginia.edu
#SBATCH --time=4-00:00:00
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6000
#SBATCH --account=csblrivanna
#SBATCH --output=/home/mac9jc/paradigm/update_universal_for_gapfilling.out

module load anaconda/5.2.0-py3.6
source activate paradigm_env
python3 denovo_reconstruction/step_3B_extend_universal.py
