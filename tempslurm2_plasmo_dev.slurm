#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=9000
#SBATCH --account=csblrivanna
#SBATCH --output=%x-%j.out

module load anaconda/5.2.0-py3.6
source activate paradigm_env

#python3 ./model_refinement/step_4_build_generic_biomass_reaction.py final_denovo_PfalciparumCD01.json
echo 'first'
#python3 ./model_refinement/step_5_orthology_based_semi_curation.py with_biomass_PfalciparumCD01.json
echo 'second'
python3 ./model_refinement/step_6_task_based_gapfilling.py ortho_PfalciparumCD01.json
echo 'third'
