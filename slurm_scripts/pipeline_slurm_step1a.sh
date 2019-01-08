#!/bin/bash
# path = "/home/mac9jc/paradigm"

# DOWNLOAD GENOMES
cd /home/mac9jc/paradigm/
Rscript ./data_acquistion/step_0_download_all_eupathDB_release41.r
echo "downloaded models"

# PREP FOR ANNOTATE GENOMES
cd /home/mac9jc/paradigm/data
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
# cd /home/mac9jc/paradigm/data/
for filename in ./*_annotated_Dec2018.fasta; do
echo "$filename"
echo "${filename:2:${#filename}-26}"
diamond blastp -d ./bigg_proteins_diamond -q $filename -o "${filename:2:${#filename}-26}_BiGG.tsv"
done
mv ./*_BiGG.tsv ./diamond_output_BiGG

# cd /home/mac9jc/paradigm/data/
for filename in ./*_annotated_Dec2018.fasta; do
echo "$filename"
echo "${filename:2:${#filename}-26}"
diamond blastp -d ./orthoMCL_proteins_diamond -q $filename -o "${filename:2:${#filename}-26}_orthoMCL.tsv"
done
mv ./*_orthoMCL.tsv ./diamond_output_orthoMCL
mv ./*_annotated_Dec2018.fasta ./genomes
echo "diamond annotation complete"

