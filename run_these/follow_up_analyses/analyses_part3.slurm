#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mac9jc@virginia.edu
#SBATCH --time=2-00:00:00
#SBATCH --partition=standard
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=6000
#SBATCH --account=csblrivanna
#SBATCH --output=%x-%j.out

module load anaconda/5.2.0-py3.6
#source activate paradigm_env
conda activate annotation_env
python3 ./model_evaluation/gene_essentiality.py
