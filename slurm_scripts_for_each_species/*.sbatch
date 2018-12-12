python3 ../denovo_reconstruction/step_3_build_de_novo_models.py '*'
python3 ../model_refinement/step_4_build_generic_biomass_reaction.py 'final_denovo_*.json'
python3 ../model_refinement/step_6_task_based_gapfilling.py 'with_biomass_denovo_*.json'
