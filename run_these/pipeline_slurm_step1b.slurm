#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mac9jc@virginia.edu
#SBATCH --time=2-00:00:00
#SBATCH --output=/home/mac9jc/paradigm/pipeline_curate_iPfal19_step1b.out
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6000
#SBATCH --account=csblrivanna


module purge
module load anaconda/5.2.0-py3.6
source activate paradigm_env
python3 ./model_refinement/step_2_curate_iPfal17.py
