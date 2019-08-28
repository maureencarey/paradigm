#!/bin/bash
# path = "/home/mac9jc/paradigm"

# DOWNLOAD GENOMES
cd ~/paradigm
echo "changed directory"
Rscript ./data_acquistion/step_0_download_all_eupathDB_release44.r
echo "downloaded genomes"

# PREP FOR ANNOTATE GENOMES
cd ~/paradigm/data
# get database sequence files to make protein database for annotating sequences
wget -O bigg_proteins.fasta 'https://github.com/cdanielmachado/carveme/raw/master/carveme/data/input/bigg_proteins.faa'
diamond makedb --in bigg_proteins.fasta -d bigg_proteins_diamond
wget -O aa_seqs_OrthoMCL_5.fasta 'http://orthomcl.org/common/downloads/release-5/aa_seqs_OrthoMCL-5.fasta'
diamond makedb --in aa_seqs_OrthoMCL_5.fasta -d orthoMCL_proteins_diamond
# get other useful files for downstream work
wget -O bigg_gprs.csv.gz 'https://github.com/cdanielmachado/carveme/raw/master/carveme/data/generated/bigg_gprs.csv.gz'
gunzip bigg_gprs.csv.gz
wget -O bigg_metabolites.txt 'http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt'
wget -O bigg_reactions.txt 'http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt'
#wget -O metanetx_chem_prop.tsv 'https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv'
#must uncomment out header of metanetx_chem_prop
echo "made annotation databases and other necessary files"

# ANNOTATE GENOMES
cd ~/paradigm/data/genomes/protein
for filename in ./*_annotatedProteins.fasta; do
echo "$filename"
echo "${filename:2:${#filename}-24}"
diamond blastp -d ~/paradigm/data/bigg_proteins_diamond -q $filename -o "${filename:2:${#filename}-26}_BiGG.tsv"
done
echo "diamond annotation against BiGG done"
mv ./*_BiGG.tsv ~/paradigm/data/diamond_output_BiGG

cd ~/paradigm/data/genomes/protein
parallel gzip ::: *.fasta
mkdir ./zipped_protein
mv *_annotatedProteins.fasta.gz ./zipped_protein
cd ~/paradigm/data/genomes/DNA
mkdir ./zipped_dna
parallel gzip ::: *.fasta
mv *.fasta.gz ./zipped_dna

# cd ~/paradigm/data/genomes/protein
# for filename in ./*_annotatedProteins.fasta; do
# diamond blastp -d ~/paradigm/data/orthoMCL_proteins_diamond -q $filename -o "${filename:2:${#filename}-26}_orthoMCL.tsv"
# done
# mv ./*_orthoMCL.tsv ~/paradigm/data/diamond_output_orthoMCL
# echo "diamond annotations complete"

