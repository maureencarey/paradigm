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
export PYTHONPATH=$EBROOTGUROBI/lib/python3.6_utf32:$PYTHONPATH
source activate paradigm_env

