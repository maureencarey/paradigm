#!/bin/bash

plasmodium_list=" PadleriG01 PbergheiANKA PbillcollinsiG01 PblacklockiG01 Pchabaudichabaudi PcoatneyiHackeri PcynomolgiB PcynomolgiM Pfalciparum3D7 Pfalciparum7G8 PfalciparumCD01 PfalciparumDd2 PfalciparumGA01 PfalciparumGB4 PfalciparumGN01 PfalciparumHB3 PfalciparumIT PfalciparumKE01 PfalciparumKH01 PfalciparumKH02 PfalciparumML01 PfalciparumSD01 PfalciparumSN01 PfalciparumTG01 PfragileNilgiri PgaboniG01 PgaboniSY75 Pgallinaceum8A PinuiSanAntonio1 PknowlesiH PknowlesiMalayanPk1A PmalariaeUG01 PovalecurtisiGH01 PpraefalciparumG01 PreichenowiCDC PreichenowiG01 PrelictumSGS1-like PvinckeipetteriCR Pvinckeivinckeivinckei Pvivax-likePvl01 PvivaxP01 PvivaxSal1 Pyoeliiyoelii17X Pyoeliiyoelii17XNL PyoeliiyoeliiYM "
falciparum=" Pfalciparum3D7 Pfalciparum7G8 PfalciparumCD01 PfalciparumDd2 PfalciparumGA01 PfalciparumGB4 PfalciparumGN01 PfalciparumHB3 PfalciparumIT PfalciparumKE01 PfalciparumKH01 PfalciparumKH02 PfalciparumML01 PfalciparumSD01 PfalciparumSN01 PfalciparumTG01 "

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
    foo=${species_string}"_stepC_for_sensitivity.sbatch"
    if [[ " $plasmodium_list " =~ .*\ $species_string\ .* ]]; then
    if [[ " $falciparum " =~ .*\ $species_string\ .* ]]; then
    cp slurm_template_lgmem.slurm $foo
    else
    cp slurm_template.slurm $foo
    fi
    echo "python3 ../model_refinement/step_6_task_based_gapfilling_for_sensitivity.py for_sensitivity_ortho_$species_string.json" >> $foo
    foo2=${species_string}"_stepD_for_sensitivity.sbatch"
    cp slurm_template.slurm $foo2
    echo "python3 ../model_refinement/step_6_task_based_gapfilling_for_sensitivity.py for_sensitivity_with_biomass_$species_string.json" >> $foo2
    else
    cp slurm_template.slurm $foo
    echo "python3 ../model_refinement/step_6_task_based_gapfilling_for_sensitivity.py for_sensitivity_with_biomass_$species_string.json" >> $foo
    fi
    sbatch $foo
    if [[ " $plasmodium_list " =~ .*\ $species_string\ .* ]]; then
    sbatch $foo2
    mv $foo2 ./slurm_scripts_for_each_species
    fi
    mv $foo ./slurm_scripts_for_each_species
done
