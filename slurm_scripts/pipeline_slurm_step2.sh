#!/bin/bash

#SBATCH --ntasks=1
#SBATCH --mail-type=ALL
#SBATCH --time=2-00:00:00
#SBATCH --partition=parallel
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=6000
#SBATCH --account=mac9jc

module load gurobi-8.0.1
module load anaconda3
source activate gurobienv
cd ./paradigm/data

# PREP LIST OF PLASMODIUM SPECIES - these models will have an extra step associated with them
plasmodium_list=" PadleriG01 PbergheiANKA PbillcollinsiG01 PblacklockiG01 Pchabaudichabaudi PcoatneyiHackeri PcynomolgiB PcynomolgiM Pfalciparum3D7 Pfalciparum7G8 PfalciparumCD01 PfalciparumDd2 PfalciparumGA01 PfalciparumGB4 PfalciparumGN01 PfalciparumHB3 PfalciparumIT PfalciparumKE01 PfalciparumKH01 PfalciparumKH02 PfalciparumML01 PfalciparumSD01 PfalciparumSN01 PfalciparumTG01 PfragileNilgiri PgaboniG01 PgaboniSY75 Pgallinaceum8A PinuiSanAntonio1 PknowlesiH PknowlesiMalayanPk1A PmalariaeUG01 PovalecurtisiGH01 PpraefalciparumG01 PreichenowiCDC PreichenowiG01 PrelictumSGS1-like PvinckeipetteriCR Pvinckeivinckeivinckei PvivaxP01 PvivaxSal1 Pyoeliiyoelii17X Pyoeliiyoelii17XNL PyoeliiyoeliiYM "

# GENERATE AND CURATE MODELS
for filename in ./diamond_output_BiGG/*_BiGG.tsv; do
    echo "$filename"
    TEXT="${filename:21}"
    species_string="${TEXT:1:${#TEXT}-10}"
    echo "$species_string"
    echo "working on: ${species_string}"
    python ../denovo_reconstruction/step_3_build_de_novo_models.py $filename
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
mv ./log_*.txt ./model_generation_logs
mkdir ./models
mv ./*.json ./models

# CHECK ALL LOGS FOR INFEASIBLE STEPS
cd .. # paradigm directory
for filename in *2018.log; do
    echo "$filename"
    python ./check_logs.py $filename
done
