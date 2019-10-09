# ParaDIGM: Parasite Database Including Genome-scale metabolic Models 

Code availability for the manuscript titled "Comparative analyses of parasites 
with a comprehensive database of genome-scale metabolic models" by Carey *et al*.

The article is submitted JOURNAL (DOI: XXX) and available on BioRxiv (see here: https://doi.org/10.1101/772467).

__Article abstract:__ Protozoan parasites cause diverse diseases with large global impacts. Research on the pathogenesis and biology of these organisms is limited by economic and experimental constraints. Accordingly, studies of one parasite are frequently extrapolated to infer knowledge about another parasite, across and within genera. Model in vitro or in vivo systems are frequently used to enhance experimental manipulability, but generally use species related to, yet distinct from, the clinically relevant causal pathogen. Characterization of functional differences among parasite species is confined to post hoc or single target studies, limiting the utility of this extrapolation approach. To address this challenge and to accelerate parasitology research broadly, we present a functional comparative analysis of 192 genomes, representing every high-quality, publicly-available protozoan parasite genome including *Plasmodium, Toxoplasma, Cryptosporidium, Entamoeba, Trypanosoma, Leishmania, Giardia,* and other species. We generated an automated metabolic network reconstruction pipeline optimized for eukaryotic organisms. These metabolic network reconstructions serve as biochemical knowledgebases for each parasite, enabling qualitative and quantitative comparisons of metabolic behavior across parasites. We identified putative differences in gene essentiality and pathway utilization to facilitate the comparison of experimental findings. This knowledgebase represents the largest collection of genome-scale metabolic models for both pathogens and eukaryotes; with this resource, we can predict species-specific functions, contextualize experimental results, and optimize selection of experimental systems for fastidious species.

##### notes

FYI some files are too large to share on github, including (but not limited to):
* genome fasta files
* compiled diamond databases
Please see .gitignore for specific file names or groups and feel free to email for access.

##### helper_functions.py
These functions are used throughout the associated scripts and called typically as 'hf'.

## data

All data in subdirectory entitled data.

##### plasmodium_orthology_conversion.csv
This is acquired from PlasmoDB.org by doing gene search> by orthology (P. falciparum 3D7) > 
add step (transform by orthology (all))

##### auxotrophies_mapping_to_genomeID.csv, auxotrophies_mapping_to_metID.csv, auxotrophies_references.xlsx
These files contain information about biochemical experiments used to gapfill each model.

##### bigg_gprs.csv
This file is part of the CarveMe package by Daniel Machado (check it out here: https://github.com/cdanielmachado/carveme)

##### Pfalciparum3D7_GeneAliases.csv
This file was obtained from PlasmoDB.org and condains the history of gene IDs used by the database.

## data acquisition

##### step_0_download_all_eupathDB_release44.r
This downloads all specified genomes from EuPathDB, release 44. The specified genomes are
listed in this script and were selected if they had fasta files for the annotated protein
sequences. Per communication with EuPathDB, the annotated protein sequences are derived 
from the official genome sequence and include a collection of manual and automatic 
ORF information. Note, there may be more genomes available in future releases.

##### step_1_genome_annotation.py
This script processes some files for visualization. The visualization is not used in this 
version of the document, but is still useful for spreadsheet-izing some of the data.

##### update_universal_reaction_set1.py
This script modifies the universal model downloaded from BiGG by calling the BiGG API to add relevant info.

## de novo reconstruction

##### step_2_build_de_novo_models.py
This builds models for all genomes from the universal model acquired from BiGG. See 
methods in publication for detail.

##### step_3B_build_de_novo_models.py
This modifies the universal model by adding all reactions from the de novo models to the universal. Why is this necessary? Well, the de novo reconstruction process generates some new reactions in the cytosol (it copies compartmentalized reactions and moves them to the cytosol). These are necessary for valid gapfilling solutions.

## model evaluation

##### gene_essentiality.py and rxn_essentiality.py
These perform essentiality simulations on a subset of models, used in Figures XXX.

##### plasmodium_history_comparison.py
This compares structural features of the gapfilled plasmodium models with and without 
orthologous transformation.

##### production_consumption.py
This is used to generate figure 5.

##### reaction_matrix_and_transporter.py
This generates a matrix of reaction and transporter presence for ungapfilled models, used in Figure 3B.

##### reaction_matrix_and_transporter_after_gf.py
This generates a matrix of reaction and transporter presence for gapfilled models, used in Figure 3E.

## model refinement

##### step_2_curate_iPfal17.py
This script curates iPfal17 to generate iPfal19.

##### step_4_build_generic_biomass_reaction.py
This builds a generic biomass reaction from curation reconstructions and adds it to all models.

##### step_5_orthology_based_semi_curation.py
This performs orthology-based semi-curation for all *Plasmodium* reconstructions using the 
curated iPfal18 reconstruction for *Plasmodium falciparum*.

##### step_6_task_based_gapfilling.py
gapfilling to generic biomass, species specific biomass (if available), and any tasks 
identified in auxotrophies_references.xlsx (if available)

##### step_7_add_annotation_obj.py
This script adds a model.annotation object with model metadata to each reconstruction. It requires the latest CobraPy.

## run these

Note: I ran the model generation pipeline on the University of Virginia's High-
Performance Computing system, Rivanna, which uses slurm job management. These were 
the scripts I ran to implement the pipeline. They are mostly self explanatory wrappers for other scripts.
Feel free to contact me with any questions.

## slurm scripts

##### check_logs.py
This is used in pipeline_cleanup.sh to parse through all model building log files to identify if any infeasible solutions were generated and writes a new log file with the problems that had infeasible solutions.

##### pipeline_slurm_step1a.sh
This is used for data acquisition in pipeline_slurm_step1.slurm. It first calls step_0_download_all_eupathDB_release44.r to download genomes and also downloads BiGG database files, and then annotates genomes (against the BiGG database and OrthoMCL) using diamond. BiGG is a project out of UCSD, see e.g. doi:
10.1093/nar/gkv1049. OrthoMCL is also a EuPathDB project, see e.g. doi: 10.1093/nar/gkj123

##### pipeline_slurm_step1b.sh
This is used for data acquisition in pipeline_slurm_step1.slurm and calls scripts in the *data acquisition* subdirectory.
This script saves annotations in table format and curates iPfal17 to generate iPfal18.

## Contact

Feedback and questions to Maureen Carey - mac9jc [at] virginia [dot] edu

## RUN ORDER ON RIVANNA

    conda create -n paradigm_env python=3.6
    conda activate paradigm_env 
    pip install -r requirements.txt
    # # manually doubel check latest EuPathDB release to see if any extra files need to run
    # # get data
    sbatch run_these/pipeline_slurm_step1.slurm
    # # curate iPfal17 model THIS CAN BE RUN CONCURRENTLY WITH PREVIOUS STEP
    sbatch run_these/pipeline_slurm_step1b.slurm
    # # add data to universal model
    sbatch run_these/update_universal_reaction_set1.slurm 
    # # make all de novo models 
    bash run_these/pipeline_auto_slurm_for_step2a.sh
    # # extend universal model by de novo models
    sbatch run_these/update_universal_reaction_set2.slurm
    # # finish making all models
    bash run_these/pipeline_auto_slurm_for_step2b.sh 
        # # TO DO: fix LmajorSD third line, remove ‘.1’, otherwise the script will fail
    # # gapfills plasmodium models prior to orthology conversion to test differences
    bash run_these/pipeline_auto_slurm_for_step2c.sh 
    # # clean things up - especially log files
    sbatch run_these/pipeline_cleanup.slurm
    # # move things to convenient locations
    mv data/*.out data/slurm_outputs
    mv data/model_modifications_* data/model_modifications
    mv data/orthology_modifications_* data/ortho_modifications
    mv data/gapfilling_additions_* data/gapfilling_additions
    mv data/percent_reactions_in_* data/percent_wrong_comp
    # # run analyses THESE CAN BE RUN CONCURRENTLY
    sbatch run_these/follow_up_analyses/analyses_part1.slurm
    sbatch run_these/follow_up_analyses/analyses_part2.slurm
    sbatch run_these/follow_up_analyses/analyses_part3.slurm
    sbatch run_these/follow_up_analyses/analyses_part4.slurm
    sbatch run_these/follow_up_analyses/analyses_part5.slurm
    sbatch run_these/follow_up_analyses/analyses_part6.slurm
    sbatch run_these/add_annotation.slurm
    # # would like to run Memote again on all models
    ## need to move all xml gf models to other directory for memote    
