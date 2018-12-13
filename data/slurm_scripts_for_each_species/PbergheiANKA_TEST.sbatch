#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --time=2-00:00:00
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6000
#SBATCH --account=csblrivanna

module load anaconda3
module load gurobi/8.0.1
source activate paradigm_env

python3 ../denovo_reconstruction/step_3_build_de_novo_models_TEST.py 'PbergheiANKA_BiGG.tsv'
python3 ../model_refinement/step_4_build_generic_biomass_reaction_TEST.py 'final_denovo_PbergheiANKA.json'
python3 ../model_refinement/step_5_orthology_based_semi_curation_TEST.py 'with_biomass_PbergheiANKA.json'
python3 ../model_refinement/step_6_task_based_gapfilling_TEST.py 'ortho_PbergheiANKA.json'
