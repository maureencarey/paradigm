#!/bin/bash

# run bioinformatic steps
cd /Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/
Rscript ./data_acquistion/step_0a_download_all_eupathDB_release41.R
echo "downloaded models"
bash ./data_acquistion/step_0b_make_diamond_database
echo "made annotation databases"
bash ./data_acquistion/step_0c_diamond_annot
echo "diamond annotation"
python ./data_acquistion/step_1_genome_annotation
echo "processed annotations"

# prep list of plasmodium species - these models will have an extra step associated with them
plasmodium_list=" PadleriG01 PbergheiANKA PbillcollinsiG01 PblacklockiG01 Pchabaudichabaudi PcoatneyiHackeri PcynomolgiB PcynomolgiM Pfalciparum3D7 PfalciparumIT PfragileNilgiri PgaboniG01 PgaboniSY75 Pgallinaceum8A PinuiSanAntonio1 PknowlesiH PknowlesiMalayanPk1A PmalariaeUG01 PovalecurtisiGH01 PpraefalciparumG01 PreichenowiCDC PreichenowiG01 PrelictumSGS1-like PvinckeipetteriCR Pvinckeivinckeivinckei PvivaxP01PvivaxSal1 Pyoeliiyoelii17X Pyoeliiyoelii17XNL PyoeliiyoeliiYM "

# finish curation of curated model
python ./model_refinement/step_3_curate_iPfal17
echo "finished iPfal17 curation"

# generate and curate models
cd /Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/data
for filename in ./diamond_output_BiGG/*_BiGG.tsv; do
    echo "$filename"
    TEXT="${filename:21}"
    species_string="${TEXT:1:${#TEXT}-10}"
    echo "$species_string"
    echo "working on: ${species_string}"
    python ../denovo_reconstruction/step_2_build_de_novo_models.py $filename
    python ../model_refinement/step_4_build_generic_biomass_reaction.py "final_denovo_$species_string.json"
    echo "finished ${species_string} de novo model generation"
    if [ "$species_string" == "$plasmodium_list" ]
    then
        python ../model_refinement/step_5_orthology_based_semi_curation.py "with_biomass_denovo_$species_string.json"
        python ../model_refinement/step_6_task_based_gapfilling.py "ortho_$species_string.json"
    else
        python ../model_refinement/step_6_task_based_gapfilling.py "with_biomass_denovo_$species_string.json"
    fi
    echo "finished reconstruction for: ${species_string}"
done
mkdir ./model_generation_logs
mv ./*.log ./model_generation_logs
mkdir ./models
mv ./*.json ./models

# check logs for infeasible steps
cd ./model_generation_logs
for filename in *2018.log; do
echo "$filename"
python ../check_logs.py $filename
done
