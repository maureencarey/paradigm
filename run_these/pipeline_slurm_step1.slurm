#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mac9jc@virginia.edu
#SBATCH --time=2-00:00:00
#SBATCH --output=/home/mac9jc/paradigm/pipeline_step1.out
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6000
#SBATCH --account=csblrivanna

#module purge
#module load gcc/7.1.0
#module load R/3.3.0
#module load diamond
#bash ./slurm_scripts/pipeline_slurm_step1a.sh

module purge
module load anaconda/5.2.0-py3.6
source activate paradigm_env
# SAVE ANNOTATIONS IN TABLE FORMAT
python3 ./data_acquistion/step_1_genome_annotation.py
