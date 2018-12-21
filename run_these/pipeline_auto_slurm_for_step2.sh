#!/bin/bash

plasmodium_list=" PadleriG01 PbergheiANKA PbillcollinsiG01 PblacklockiG01 Pchabaudichabaudi PcoatneyiHackeri PcynomolgiB PcynomolgiM Pfalciparum3D7 Pfalciparum7G8 PfalciparumCD01 PfalciparumDd2 PfalciparumGA01 PfalciparumGB4 PfalciparumGN01 PfalciparumHB3 PfalciparumIT PfalciparumKE01 PfalciparumKH01 PfalciparumKH02 PfalciparumML01 PfalciparumSD01 PfalciparumSN01 PfalciparumTG01 PfragileNilgiri PgaboniG01 PgaboniSY75 Pgallinaceum8A PinuiSanAntonio1 PknowlesiH PknowlesiMalayanPk1A PmalariaeUG01 PovalecurtisiGH01 PpraefalciparumG01 PreichenowiCDC PreichenowiG01 PrelictumSGS1-like PvinckeipetteriCR Pvinckeivinckeivinckei PvivaxP01 PvivaxSal1 Pyoeliiyoelii17X Pyoeliiyoelii17XNL PyoeliiyoeliiYM "

cd /home/mac9jc/paradigm/data
mkdir ./slurm_scripts_for_each_species

for filename in ./diamond_output_BiGG/*_BiGG.tsv; do
    echo "$filename"
    file_without_pre="${filename##*/}"
    echo "$file_without_pre"
    file_without_ext="${file_without_pre%.*}"
    echo "$file_without_ext"
    species_string="${file_without_ext:0:${#file_without_ext}-5}"
    echo "$species_string"
    #construct name for sbatch file being generated
    foo=${species_string}".sbatch"
    cp slurm_template.slurm $foo
    echo "python3 ../denovo_reconstruction/step_3_build_de_novo_models.py $file_without_pre" >> $foo
    echo "python3 ../model_refinement/step_4_build_generic_biomass_reaction.py final_denovo_$species_string.json" >> $foo
    if [[ " $plasmodium_list " =~ .*\ $species_string\ .* ]]; then
    echo "python3 ../model_refinement/step_5_orthology_based_semi_curation.py with_biomass_$species_string.json" >> $foo
    echo "python3 ../model_refinement/step_6_task_based_gapfilling.py ortho_$species_string.json" >> $foo
    else
    echo "python3 ../model_refinement/step_6_task_based_gapfilling.py with_biomass_$species_string.json" >> $foo
    fi
    sbatch $foo
    mv $foo ./slurm_scripts_for_each_species
done

