import cobra
import os
import pandas as pd
import numpy as np
import glob
from cobra.core import Gene, Metabolite, Reaction

from datetime import datetime

day = datetime.now().strftime('%d_%m_%Y')

## ##### get essential
os.chdir("/home/mac9jc/paradigm/models")
essentiality_screen_models = dict()
essentiality_screen_models['TgondiiGT1'] = cobra.io.read_sbml_model('gf_no_ortho_TgondiiGT1.xml')
essentiality_screen_models['TgondiiME49'] = cobra.io.read_sbml_model('gf_no_ortho_TgondiiME49.xml')
essentiality_screen_models['Pfalciparum3D7'] = cobra.io.read_sbml_model('gf_Pfalciparum3D7.xml')
essentiality_screen_models['PfalciparumDd2'] = cobra.io.read_sbml_model('gf_PfalciparumDd2.xml')
essentiality_screen_models['PbergheiANKA'] = cobra.io.read_sbml_model('gf_PbergheiANKA.xml')

essentiality_screen_models['TgondiiGT1_for_sensitivity'] = cobra.io.read_sbml_model('for_sensitivity_gf_no_ortho_TgondiiGT1.xml')
essentiality_screen_models['TgondiiME49_for_sensitivity'] = cobra.io.read_sbml_model('for_sensitivity_gf_no_ortho_TgondiiME49.xml')
essentiality_screen_models['Pfalciparum3D7_for_sensitivity'] = cobra.io.read_sbml_model('for_sensitivity_gf_Pfalciparum3D7.xml')
essentiality_screen_models['PfalciparumDd2_for_sensitivity'] = cobra.io.read_sbml_model('for_sensitivity_gf_PfalciparumDd2.xml')
essentiality_screen_models['PbergheiANKA_for_sensitivity'] = cobra.io.read_sbml_model('for_sensitivity_gf_PbergheiANKA.xml')

gene_essentiality_screen_results_raw= dict()
gene_essentiality_screen_results_interpreted = dict()

for species, model in essentiality_screen_models.items():
    raw_results = dict()
    interpreted_results = dict()

    print(species+', gene essenitality screen')
    
    use_second_biomass = False
    if species in ['TgondiiGT1','TgondiiME49','TgondiiGT1_for_sensitivity','TgondiiME49_for_sensitivity']:
        model.objective = 'generic_biomass'
        model.reactions.get_by_id('generic_biomass').upper_bound = 1000.
        model.reactions.get_by_id('generic_biomass').lower_bound = 0.  
    else:
        model.objective = "generic_biomass"
        use_second_biomass = True
        model.reactions.get_by_id('generic_biomass').upper_bound = 1000.
        model.reactions.get_by_id('generic_biomass').lower_bound = 0.
        # don't accidentally use other biomass reaction
        model.reactions.get_by_id('biomass').upper_bound = 0.
        model.reactions.get_by_id('biomass').lower_bound = 0.

    max_biomass = model.slim_optimize()
    print(max_biomass)
    if max_biomass < 0.01:
        print('model doesnt grow')
        print(species)
          
    # knockout and record growth
    for gene in model.genes:
        with model as cobra_model:
            cobra.manipulation.delete_model_genes(cobra_model, [gene.id], cumulative_deletions=True)
            f = cobra_model.slim_optimize()
            if max_biomass < 0.01:
       	       	print(species)
            else:
                if f < 0.1*max_biomass:
                    interpreted_results[gene.id] = 'lethal'
                else:
                    interpreted_results[gene.id] = 'nonlethal'
                raw_results[gene.id] = f/max_biomass
            cobra.manipulation.undelete_model_genes(cobra_model)

    # save
    if not use_second_biomass:
        print('only one biomass reaction')
        gene_essentiality_screen_results_raw[species] = raw_results
        gene_essentiality_screen_results_interpreted[species] = interpreted_results
        
        os.chdir("/home/mac9jc/paradigm/data")
        pd.DataFrame.from_dict(gene_essentiality_screen_results_raw[species], orient='index').to_csv("/home/mac9jc/paradigm/data/results/gene_essentiality/for_sensitivity_gene_essentiality_matrix_{}_{}.csv".format(species,day))
        pd.DataFrame.from_dict(gene_essentiality_screen_results_interpreted[species], orient='index').to_csv("/home/mac9jc/paradigm/data/results/gene_essentiality/for_sensitivity_gene_essentiality_matrix_interpreted_{}_{}.csv".format(species,day))

    else:
       	print('wroking on second biomass reactions')
        gene_essentiality_screen_results_raw[species+'_generic_biomass'] = raw_results
        gene_essentiality_screen_results_interpreted[species+'_generic_biomass'] = interpreted_results
        
        pd.DataFrame.from_dict(gene_essentiality_screen_results_raw[species+'_generic_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/results/gene_essentiality/for_sensitivity_gene_essentiality_matrix_{}_{}.csv".format(species+'_generic_biomass',day))
        pd.DataFrame.from_dict(gene_essentiality_screen_results_interpreted[species+'_generic_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/results/gene_essentiality/for_sensitivity_gene_essentiality_matrix_interpreted_{}_{}.csv".format(species+'_generic_biomass',day))

        print(species+', gene essenitality screen')
        # set objective
        model.objective = "biomass"

        # don't accidentally use other biomass reaction
        model.reactions.get_by_id('biomass').upper_bound = 1000.
        model.reactions.get_by_id('biomass').lower_bound = 0.
        model.reactions.get_by_id('generic_biomass').upper_bound = 0.
        model.reactions.get_by_id('generic_biomass').lower_bound = 0.

        max_biomass = model.slim_optimize()
        if max_biomass < 0.01:
            print('model doesnt grow')
            print(species)

        # knockout and record growth
        for gene in model.genes:
            with model as cobra_model:
                cobra.manipulation.delete_model_genes(cobra_model, [gene.id], cumulative_deletions=True)
                f = cobra_model.slim_optimize()
                if max_biomass < 0.01:
       	       	    print(species)
                else:
                    if f < 0.1*max_biomass:
                        interpreted_results[gene.id] = 'lethal'
                    else:
                        interpreted_results[gene.id] = 'nonlethal'
                    raw_results[gene.id] = f/max_biomass
                cobra.manipulation.undelete_model_genes(cobra_model)

        # save
        gene_essentiality_screen_results_raw[species+'_species_biomass'] = raw_results
        gene_essentiality_screen_results_interpreted[species+'_species_biomass'] = interpreted_results

        pd.DataFrame.from_dict(gene_essentiality_screen_results_raw[species+'_species_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/results/gene_essentiality/for_sensitivity_gene_essentiality_matrix_{}_{}.csv".format(species+'_species_biomass',day))
        pd.DataFrame.from_dict(gene_essentiality_screen_results_interpreted[species+'_species_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/results/gene_essentiality/for_sensitivity_gene_essentiality_matrix_interpreted_{}_{}.csv".format(species+'_species_biomass',day))



