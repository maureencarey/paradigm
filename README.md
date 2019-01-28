# models for eukaryotic organisms and semicuration approach

Code availability for the manuscript titled "title" by Carey *et al*.

The article is submitted JOURNAL (DOI: XXX) and available on BioRxiv (see here: XXX).

FYI some files are to large to share on github, including (but not limited to):
* genome fasta files
* compiled diamond databases
Please see .gitignore for specific file names or groups and feel free to email for access.

## data

All data in subdirectory entitled data.

## data acquisition

### step_0_download_all_eupathDB_release41.r

This downloads all specified genomes from EuPathDB, release 41. The specified genomes are
listed in this script and were selected if they had fasta files for the annotated protein
sequences. Per communication with EuPathDB, the annotated protein sequences are derived 
from the official genome sequence and include a collection of manual and automatic 
ORF information. Note, there may be more genomes available in future releases.

### step_1_genome_annotation.py

This script processes some files for visualization. The visualization is not used in this version of the document, but is still useful for spreadsheet-izing some of the data.

## de novo reconstruction

### helper_functions_1.py

This functions are used throughout the associated scripts and called typically as 'hf'.


### step_2_build_de_novo_models.py

This builds models for all genomes from the universal model acquired from BiGG. See methods in publication for detail.

## model evaluation

### gene_essentiality.py and rxn_essentiality.py

These perform essentiality simulations on a subset of models, used in Figures XXX.

### plasmodium_history_comparison.py

This compares structural features of the gapfilled plasmodium models with and without orthologous transformation.

### production_consumption.py

This is used to generate figure XXX.

### reaction_matrix_and_transporter.py

This generates a matrix of reaction and transporter presence for ungapfilled models, used in Figure XXX.

## model refinement

### helper_functions_2.py

This functions are used throughout the associated scripts and called typically as 'hf'.

### step_4_build_generic_biomass_reaction.py

This builds a generic biomass reaction from curation reconstrucitons and adds it to all models.


### step_5_orthology_based_semi_curation.py

This performs orthology-based semi-curation for all *Plasmodium* reconstructions using the curated iPfal18 reconstruction
for *Plasmodium falciparum*.

### step_6_task_based_gapfilling.py

gapfilling to generic biomass, species specific biomass (if available), and any tasks idenified in FILENAME (if available)


### step_6_task_based_gapfilling_plasmodium.py

This is the same as step_6_task_based_gapfilling.py but save the files with obviously different names when gapfilling the plasmodium models WITHOUT orthologous transformation.

### TO DO: remove futile cycles

## models

Interestingly enough, this directory contains the models.


## run these

Note: I ran the model generation pipeline on the University of Virginia's High-Performance Computing system, Rivanna, which uses Slurm job management. These were the scripts I ran to implement the pipeline. They are mostly 'wrappers' for other scripts.

### follow up analyses

This directory contains the analyses scripts - all are self explanatory and point towards the python scripts contained in other subdirectories.

### pipeline_auto_slurm_for_plasmodium.sh

This file generates the sbatch scripts for gapfilling the de novo plasmodium reconstructions WITHOUT the orthology transformation.

### pipeline_auto_slurm_for_step2.sh

This file generates the sbatch scripts for generating all reconstructions. Only the plasmodium recosntructions use orthology transformation.

### pipeline_cleanup.sh

This script is used to interpret log files from the reconstruction process.

### pipeline_slurm_step1.slurm

This slurm script runs the data acquisition scripts in subdirecotry *slurm_scripts*.

### update_universal_reaction_set.py

This script removes the biomass reactions from the universal reaction set from BiGG and adds some annotation info that is stored on BiGG to the metabolites and reactions.

## slurm scripts

### check_logs.py

This is used in pipeline_cleanup.sh to parse through all model building log files to identify if any infeasible solutions were generated and writes a new log file with the problems that had infeasible solutions.

### pipeline_slurm_step1a.sh

This is used for data acquisition in pipeline_slurm_step1.slurm

It first calls step_0_download_all_eupathDB_release41.r to download genomes and also downloads BiGG database files, and then annotates genomes (against the BiGG database and OrthoMCL) using diamond. BiGG is a project out of UCSD, see e.g. doi:
10.1093/nar/gkv1049. OrthoMCL is also a EuPathDB project, see e.g. doi: 10.1093/nar/gkj123

### pipeline_slurm_step1b.sh

This is used for data acquisition in pipeline_slurm_step1.slurm and calls scripts in the *data acquisition* subdirectory.

This script saves annotations in table format and curates iPfal17 to generate iPfal18.


## NOT DONE

### metadata_conversion.py


## Contact

Feedback and questions to Maureen Carey - mac9jc [at] virginia [dot] edu

## RUN ORDER ON RIVANNA

    # # manually doubel check latest EuPathDB release to see if any extra files need to run
    # # get data
    sbatch ./run_these/pipeline_slurm_step1.slurm
    module load anaconda/5.2.0-py3.6
    # # would like to run Memote here on iPfal18
    python 3 ./run_these/update_universal_reaction_set.py
    # # would like to run Memote again on iPfal18
    # # make all models
    bash ./run_these/pipeline_auto_slurm_for_step2.sh
    # # TO DO: fix LmajorSD third line, remove ‘.1’, otherwise the script will fail
    # # clean things up - especially log files
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
    # # would like to run Memote again on all models
    
