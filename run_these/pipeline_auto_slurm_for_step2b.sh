#!/bin/bash

plasmodium_list=" PadleriG01 PbergheiANKA PbillcollinsiG01 PblacklockiG01 Pchabaudichabaudi PcoatneyiHackeri PcynomolgiB PcynomolgiM Pfalciparum3D7 Pfalciparum7G8 PfalciparumCD01 PfalciparumDd2 PfalciparumGA01 PfalciparumGB4 PfalciparumGN01 PfalciparumHB3 PfalciparumIT PfalciparumKE01 PfalciparumKH01 PfalciparumKH02 PfalciparumML01 PfalciparumSD01 PfalciparumSN01 PfalciparumTG01 PfragileNilgiri PgaboniG01 PgaboniSY75 Pgallinaceum8A PinuiSanAntonio1 PknowlesiH PknowlesiMalayanPk1A PmalariaeUG01 PovalecurtisiGH01 PpraefalciparumG01 PreichenowiCDC PreichenowiG01 PrelictumSGS1-like PvinckeipetteriCR Pvinckeivinckeivinckei PvivaxP01 PvivaxSal1 Pyoeliiyoelii17X Pyoeliiyoelii17XNL PyoeliiyoeliiYM "

cd /home/mac9jc/paradigm/data
# mkdir ./slurm_scripts_for_each_species

for filename in ./diamond_output_BiGG/*_BiGG.tsv; do
    echo "$filename"
    file_without_pre="${filename##*/}"
    echo "$file_without_pre"
    file_without_ext="${file_without_pre%.*}"
    echo "$file_without_ext"
    species_string="${file_without_ext:0:${#file_without_ext}-5}"
    echo "$species_string"
    #construct name for sbatch file being generated
    foo=${species_string}"_stepB.sbatch"
    cp slurm_template.slurm $foo
    if [[ " $plasmodium_list " =~ .*\ $species_string\ .* ]]; then
    echo "python3 ../model_refinement/step_5_orthology_based_semi_curation.py with_biomass_$species_string.json" >> $foo
    echo "python3 ../model_refinement/step_6_task_based_gapfilling.py ortho_$species_string.json" >> $foo
    else
    echo "python3 ../model_refinement/step_6_task_based_gapfilling.py with_biomass_$species_string.json" >> $foo
    fi
    sbatch $foo
    mv $foo ./slurm_scripts_for_each_species
done

cd /home/mac9jc/paradigm/data
mkdir ./slurm_scripts_for_without_ortho_plasmodium

for filename in with_biomass_PadleriG01.json with_biomass_PbergheiANKA.json with_biomass_PbillcollinsiG01.json with_biomass_PblacklockiG01.json with_biomass_Pchabaudichabaudi.json with_biomass_PcoatneyiHackeri.json with_biomass_PcynomolgiB.json with_biomass_PcynomolgiM.json with_biomass_Pfalciparum3D7.json with_biomass_Pfalciparum7G8.json with_biomass_PfalciparumCD01.json with_biomass_PfalciparumDd2.json with_biomass_PfalciparumGA01.json with_biomass_PfalciparumGB4.json with_biomass_PfalciparumGN01.json with_biomass_PfalciparumHB3.json with_biomass_PfalciparumIT.json with_biomass_PfalciparumKE01.json with_biomass_PfalciparumKH01.json with_biomass_PfalciparumKH02.json with_biomass_PfalciparumML01.json with_biomass_PfalciparumSD01.json with_biomass_PfalciparumSN01.json with_biomass_PfalciparumTG01.json with_biomass_PfragileNilgiri.json with_biomass_PgaboniG01.json with_biomass_PgaboniSY75.json with_biomass_Pgallinaceum8A.json with_biomass_PinuiSanAntonio1.json with_biomass_PknowlesiH.json with_biomass_PknowlesiMalayanPk1A.json with_biomass_PmalariaeUG01.json with_biomass_PovalecurtisiGH01.json with_biomass_PpraefalciparumG01.json with_biomass_PreichenowiCDC.json with_biomass_PreichenowiG01.json with_biomass_PrelictumSGS1-like.json with_biomass_PvinckeipetteriCR.json with_biomass_Pvinckeivinckeivinckei.json with_biomass_PvivaxP01.json with_biomass_PvivaxSal1.json with_biomass_Pyoeliiyoelii17XNL.json with_biomass_Pyoeliiyoelii17X.json with_biomass_PyoeliiyoeliiYM.json; do
    echo "$filename"
    foo=${filename}".sbatch"
    cp slurm_template_plasmodium.slurm $foo
    echo "python3 ../model_refinement/step_6_task_based_gapfilling_plasmodium.py /home/mac9jc/paradigm/models/$filename" >> $foo
    sbatch $foo
    mv $foo ./slurm_scripts_for_without_ortho_plasmodium
done



