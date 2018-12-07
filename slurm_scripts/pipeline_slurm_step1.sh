#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --time=2-00:00:00
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6000
#SBATCH --account=MY_ACCOUNT

module load gurobi-7.0.2
module load anaconda3
source activate gurobienv
./paradigm

# run bioinformatic steps
R ./data_acquistion/step_0a_download_all_eupathDB_release.R # update to 40
echo "downloaded models"
bash ./data_acquistion/step_0b_make_diamond_database
echo "made annotation databases"
bash ./data_acquistion/step_0c_diamond_annot
echo "diamond annotation"
python ./data_acquistion/step_1_genome_annotation
echo "processed annotations"

# finish curation of curated model
python3 ./model_refinement/step_3_curate_iPfal17
echo "finished iPfal17 curation"

