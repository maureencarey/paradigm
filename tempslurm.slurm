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

python3 ../denovo_reconstruction/step_3_build_de_novo_models.py TgondiiME49_BiGG.tsv
#python3 ./model_refinement/step_4_build_generic_biomass_reaction.py final_denovo_TgondiiME49.json
#python3 ./model_refinement/step_6_task_based_gapfilling.py with_biomass_TgondiiME49.json
