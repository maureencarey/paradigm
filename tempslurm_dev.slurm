#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --partition=dev
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=9000
#SBATCH --account=csblrivanna
#SBATCH --output=%x-%j.out

module load anaconda/5.2.0-py3.6
source activate paradigm_env

#python3 ./denovo_reconstruction/step_3_build_de_novo_models.py AalgeraePRA109_BiGG.tsv
#echo 'first script'
python3 ./model_refinement/step_4_build_generic_biomass_reaction.py final_denovo_AalgeraePRA109.json
echo 'second'
python3 ./model_refinement/step_6_task_based_gapfilling.py with_biomass_AalgeraePRA109.json
echo 'third script'
