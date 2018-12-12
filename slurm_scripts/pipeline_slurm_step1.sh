#!/bin/bash
# path = "/home/mac9jc/"
cd ./paradigm

# DOWNLOAD GENOMES
module load R/3.4.3
Rscript ./data_acquistion/step_0_download_all_eupathDB_release41.R
echo "downloaded models"

# PREP FOR ANNOTATE GENOMES
cd ./paradigm/data
module load diamond
# get database sequence files to make protein database for annotating sequences
wget -O bigg_proteins.fasta 'https://github.com/cdanielmachado/carveme/raw/master/carveme/data/input/bigg_proteins.faa'
diamond makedb --in bigg_proteins.fasta -d bigg_proteins_diamond
wget -O aa_seqs_OrthoMCL_5.fasta 'http://orthomcl.org/common/downloads/release-5/aa_seqs_OrthoMCL-5.fasta'
diamond makedb --in aa_seqs_OrthoMCL_5.fasta -d orthoMCL_proteins_diamond
# get other useful files for downstream work
wget -O bigg_gprs.csv.gz 'https://github.com/cdanielmachado/carveme/raw/master/carveme/data/generated/bigg_gprs.csv.gz'
gunzip bigg_gprs.csv.gz
wget -O bigg_metabolites.txt 'http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt'
echo "made annotation databases and other necessary files"

# ANNOTATE GENOMES
for filename in ./*_annotated_Dec2018.fasta; do
echo "$filename"
echo "${filename:2:${#filename}-26}"
diamond blastp -d ./bigg_proteins_diamond -q $filename -o "${filename:2:${#filename}-26}_BiGG.tsv"
done
mkdir ./diamond_output_BiGG
mv ./*_BiGG.tsv ./diamond_output_BiGG

for filename in ./*_annotated_Oct2018.fasta; do
echo "$filename"
echo "${filename:2:${#filename}-26}"
diamond blastp -d ./orthoMCL_proteins_diamond -q $filename -o "${filename:2:${#filename}-26}_orthoMCL.tsv"
done
mkdir ./diamond_output_orthoMCL
mv ./*_orthoMCL.tsv ./diamond_output_orthoMCL
mkdir ./genomes
mv ./*_annotated_Dec2018.fasta ./genomes
echo "diamond annotation complete"

# SAVE ANNOTATIONS IN TABLE FORMAT
cd .. # paradigm directory
module load anaconda3
source activate paradigm_env
python3 ./data_acquistion/step_1_genome_annotation
echo "processed annotations"

# FINISH CURATION OF iPfal17
module load gurobi/8.0.1
python3 ./model_refinement/step_2_curate_iPfal17
echo "finished iPfal17 curation"

