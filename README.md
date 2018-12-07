# models for eukaryotic organisms and semicuration approach

Code availability for the manuscript titled "title" by Carey *et al*.

The article is submitted JOURNAL (DOI: XXX) and available on BioRxiv (see here: XXX).

Reminder, might need to do: chmod +x /path/to/yourscript.sh to execute some bash scripts.
## data

All data in subdirectory entitled data.

## data acquisition

### download_all_eupathDB_release39

This downloads all specified genomes from EuPathDB, release 39. The specified genomes are 
listed in this script and were selected if they had fasta files for the annotated protein
sequences. Per communication with EuPathDB, the annotated protein sequences are derived 
from the official genome sequence and include a collection of manual and automatic 
ORF information. Note, there may be more genomes available in future releases. 
outputs:

### genome_annotation (using diamond_annot and make_diamond_database)

This annotates amino acid sequences obtained from download_all_eupathDB_release39 using 
two databases, BiGG and OrthoMCL. BiGG is a project out of UCSD, see e.g. doi:
10.1093/nar/gkv1049. OrthoMCL is also a EuPathDB project, see e.g. doi: 10.1093/nar/gkj123
inputs:
outputs:

### make_diamond_database.sh

This generates a database from fasta files. This output is necessary for genome_annotation
inputs:
outputs:

## dependencies
Required licenses and files for implementation.
Includes:

## helper functions

### helper_functions.py

This functions are used throughout the associated scripts and called typically as 'hf'.

### extract_from_ensemble.py

This script will extract an individual model from an ensemble of models.
Requirements: Medusa (https://github.com/gregmedlock/Medusa). *Alternatively, you can
access individual models in the following subdirectory: data/models/. Models are grouped
by EuPathDB database.*

### metadata_conversion.py

NOT DONE

## model evaluation

NOT DONE

### gene_essentiality.ipnb

## model generation

### build_de_novo_models.py

This builds models for all genomes acquired from download_all_eupathDB_release39.R from
the universal model acquired from BiGG.

Input:
   - annotation files from genome_annotation.py
   - universal_model (acquired here: )
Output:
   - ungapfilled models for all genomes

### build_generic_biomass_reaction.py

This builds a generic biomass reaction from curation reconstrucitons and adds it to all models.

## model refinement

### gapfilling.py

### preliminary_iPfal17_curation.ipynb CURRENTLY NAMED iPfal17_curation_part1.ipyb

This implements a few curation steps (see Supplemental Table XXX) to the Plasmodium
falciparum model published in Carey et al. BMC Genomics 2017.
Input:
iPfal17 (filename, model from Carey et al. BMC Genomics 2017)
curation file (filename, list of curation in Supplemental Table XXX)
Output:
iPfal18

### orthology_based_semi_curation.py

This performs orthology-based semi-curation for all *Plasmodium* reconstructions using the curated iPfal18 reconstruction
for *Plasmodium falciparum*.

## other analyses

## Scripts for data acquisition

### download_all_eupathDB_release39.R

This script will download all genomes on EuPathDB and was used to download genomes from 
release 39. Requirements: tidyverse, RCurl, and a whole lot of memory (about XXX gb).

## rmarkdown

This includes files to compile the manuscript.

## visualization

### upsetR_models.R

This script generates figure XXX. Requirements: UpSetR.

### transporter_similarity.R

This script generates figure XXX. Requirements:

## RUN ORDER

0. DOwnload FILENAMES
1. download_all_eupathDB_release39.R
2. make_diamond_database.sh
3.


## Contact

Feedback and questions to Maureen Carey - mac9jc [at] virginia [dot] edu
