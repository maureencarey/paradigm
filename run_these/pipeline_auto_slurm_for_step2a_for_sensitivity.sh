#!/bin/bash

cd /home/mac9jc/paradigm/data
# mkdir ./slurm_scripts_for_each_species

FILES="diamond_output_BiGG/PbergheiANKA_BiGG.tsv
diamond_output_BiGG/Pfalciparum3D7_BiGG.tsv
diamond_output_BiGG/PfalciparumDd2_BiGG.tsv
diamond_output_BiGG/TgondiiME49_BiGG.tsv
diamond_output_BiGG/TgondiiGT1_BiGG.tsv"

for filename in $FILES; do
    echo "$filename"
    file_without_pre="${filename##*/}"
    echo "$file_without_pre"
    file_without_ext="${file_without_pre%.*}"
    echo "$file_without_ext"
    species_string="${file_without_ext:0:${#file_without_ext}-5}"
    echo "$species_string"
    #construct name for sbatch file being generated
    foo=${species_string}"_stepA_for_sensitivity.sbatch"
    cp slurm_template.slurm $foo
    echo "python3 /home/mac9jc/paradigm/denovo_reconstruction/step_3_build_de_novo_models_for_sensitivity.py $file_without_pre" >> $foo
    echo "python3 /home/mac9jc/paradigm//model_refinement/step_4_build_generic_biomass_reaction_for_sensitivity.py for_sensitivity_denovo_$species_string.json" >> $foo
    sbatch $foo
    mv $foo ./slurm_scripts_for_each_species
done

