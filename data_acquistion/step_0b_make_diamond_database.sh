cd /Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/data

# get database sequence files to make protein database for annotating sequences
wget -O bigg_proteins.fasta 'https://github.com/cdanielmachado/carveme/raw/master/carveme/data/input/bigg_proteins.faa'
diamond makedb --in bigg_proteins.fasta -d bigg_proteins_updated
wget -O aa_seqs_OrthoMCL_5.fasta 'http://orthomcl.org/common/downloads/release-5/aa_seqs_OrthoMCL-5.fasta'
diamond makedb --in aa_seqs_OrthoMCL_5.fasta -d orthoMCL_proteins

# get other useful files for downstream work
wget -O bigg_gprs.csv.gz 'https://github.com/cdanielmachado/carveme/raw/master/carveme/data/generated/bigg_gprs.csv.gz'
gunzip bigg_gprs.csv.gz
wget -O bigg_metabolites.txt 'http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt'
