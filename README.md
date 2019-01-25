# models for eukaryotic organisms and semicuration approach

Code availability for the manuscript titled "title" by Carey *et al*.

The article is submitted JOURNAL (DOI: XXX) and available on BioRxiv (see here: XXX).

Reminder, might need to do: chmod +x /path/to/yourscript.sh to execute some bash scripts.

FYI some files are to large to share on github, including (but not limited to):
    - genome fasta files
    - compiled diamond databases
please see .gitignore for specific file names or groups and feel free to email for access.

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

## de novo reconstruction

### helper_functions_1.py

This functions are used throughout the associated scripts and called typically as 'hf'.

### step_2_build_de_novo_models.py

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

### helper_functions_2.py

This functions are used throughout the associated scripts and called typically as 'hf'.

### step_3_curate_iPfal17.py

This implements a few curation steps (see Supplemental Table XXX) to the Plasmodium
falciparum model published in Carey et al. BMC Genomics 2017.
Input:
iPfal17 (filename, model from Carey et al. BMC Genomics 2017)
curation file (filename, list of curation in Supplemental Table XXX)
Output:
iPfal18

### step_4_build_generic_biomass_reaction.py

### step_5_orthology_based_semi_curation.py

This performs orthology-based semi-curation for all *Plasmodium* reconstructions using the curated iPfal18 reconstruction
for *Plasmodium falciparum*.

### step_6_task_based_gapfilling.py

gapfilling to generic biomass, species specific biomass (if available), and any tasks idenified in FILENAME (if available)


## other analyses

## Scripts for data acquisition

### download_all_eupathDB_release39.R

This script will download all genomes on EuPathDB and was used to download genomes from 
release 39. Requirements: tidyverse, RCurl, and a whole lot of memory (about XXX gb).

## rmarkdown

This includes files to compile the manuscript.

## Slurm scripts

Note: I ran the model generation pipeline on the University of Virginia's High-Performance Computing system, Rivanna, which uses Slurm job management. If run on a different system, I recommend using pipeline.sh. The functionality is identical in the group of slurm files and in pipeline.sh.

### pipeline.sh

This compiles all the information presented in the next three scripts for running all model building scripts.

### pipeline_slurm_step1.slurm
### pipeline_slurm_step2.slurm
### pipeline_slurm_step3.slurm

### check_logs.sh

This file parse through a model building log files to identify if any infeasible solutions were generated and writes a new log file with the problems that had infeasible solutions.


## Contact

Feedback and questions to Maureen Carey - mac9jc [at] virginia [dot] edu

## RUN ORDER ON RIVANNA

    # # manually doubel check latest EuPathDB release to see if any extra files need to run
    # # get data
    sbatch ./run_these/pipeline_slurm_step1.slurm
    # # make all models
    bash ./run_these/pipeline_auto_slurm_for_step2.sh
    # # TO DO: fix LmajorSD third line, remove ‘.1’, otherwise the script will fail
    # # clean things up - especially log files
    module load anaconda/5.2.0-py3.6
    bash ./run_these/pipeline_cleanup.sh
	# infeasible
		# step6_EhistolyticaHM1IMSS_07_01_2019.log
		# step6_EhistolyticaHM3IMSS_07_01_2019.log
		# step6_EhistolyticaHM1IMSS-B_07_01_2019.log
		# step6_EhistolyticaKU27_07_01_2019.log
		# step6_LmexicanaMHOMGT2001U1103_07_01_2019.log
		# step6_PbergheiANKA_07_01_2019.log
		# step6_TbruceigambienseDAL972_07_01_2019.log
		# step6_TbruceiLister427_07_01_2019.log
		# step6_TbruceiTREU927_07_01_2019.log
		# step6_TgondiiARI_07_01_2019
		# step6_TgondiiFOU_07_01_2019
		# step6_TgondiiGAB2-2007-GAL-DOM2_07_01_2019
		# step6_TgondiiGT1_07_01_2019
		# step6_TgondiiMAS_07_01_2019
		# step6_TgondiiME49_07_01_2019
		# step6_TgondiiVAND_07_01_2019
			# pFBA gapfilling for DM_mal__D_c is infeasible!
    # # gapfill plasmodium models prior to orthology conversion to test differences
    bash ./run_these/pipeline_auto_slurm_for_plasmodium.sh
    # # move things to convenient locations
    # mkdir ./slurm_outputs
    mv *.out ./slurm_outputs
    # mkdir ./model_modifications
    mv model_modifications_* ./model_modifications
    # mkdir ./ortho_modifications
    mv orthology_modifications_* ./ortho_modifications
    # mkdir ./gapfilling_additions
    mv gapfilling_additions_* ./gapfilling_additions
    # mkdir ./percent_wrong_comp
    mv percent_reactions_in_* ./percent_wrong_com
    # # run analyses
    sbatch ./run_these/follow_up_analyses/analyses_part1.slurm
    sbatch ./run_these/follow_up_analyses/analyses_part2.slurm
    sbatch ./run_these/follow_up_analyses/analyses_part3.slurm
    sbatch ./run_these/follow_up_analyses/analyses_part4.slurm
    sbatch ./run_these/follow_up_analyses/analyses_part5.slurm
    
