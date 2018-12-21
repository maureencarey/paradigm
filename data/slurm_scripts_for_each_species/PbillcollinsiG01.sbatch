#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --time=2-00:00:00
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6000
#SBATCH --account=csblrivanna

module load anaconda/5.2.0-py3.6
module load gurobi/8.0.1
source activate paradigm_env

python3 ../denovo_reconstruction/step_3_build_de_novo_models.py PbillcollinsiG01_BiGG.tsv
python3 ../model_refinement/step_4_build_generic_biomass_reaction.py 'final_denovo_PbillcollinsiG01.json'
python3 ../model_refinement/step_5_orthology_based_semi_curation.py 'with_biomass_PbillcollinsiG01.json'
python3 ../model_refinement/step_6_task_based_gapfilling.py 'ortho_PbillcollinsiG01.json'
