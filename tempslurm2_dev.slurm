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

#python3 ./model_refinement/step_4_build_generic_biomass_reaction.py final_denovo_CparvumIowaII.json
#echo 'second'
python3 ./model_refinement/step_6_task_based_gapfilling.py with_biomass_CparvumIowaII.json
#echo 'third script'
