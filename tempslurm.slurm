#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --time=2-00:00:00
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=9000
#SBATCH --account=csblrivanna

module load anaconda/5.2.0-py3.6
source activate paradigm_env

python3 ../denovo_reconstruction/step_3_build_de_novo_models.py AalgeraePRA109_BiGG.tsv
