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
essentiality_screen_models['PbergheiANKA'] = cobra.io.read_sbml_model('gf_PbergheiANKA.xml')
essentiality_screen_models['PcynomolgiB'] = cobra.io.read_sbml_model('gf_PcynomolgiB.xml')
essentiality_screen_models['PvivaxSal1'] = cobra.io.read_sbml_model('gf_PvivaxSal1.xml')
essentiality_screen_models['ChominisTU502_2012'] = cobra.io.read_sbml_model('gf_no_ortho_ChominisTU502_2012.xml')
essentiality_screen_models['CparvumIowaII'] = cobra.io.read_sbml_model('gf_no_ortho_CparvumIowaII.xml')
essentiality_screen_models['PknowlesiH'] = cobra.io.read_sbml_model('gf_PknowlesiH.xml')
essentiality_screen_models['PcynomolgiB'] = cobra.io.read_sbml_model('gf_PcynomolgiB.xml')

essentiality_screen_models['Pfalciparum3D7_without_ortho'] = cobra.io.read_sbml_model('gf_no_ortho_Pfalciparum3D7.xml')
essentiality_screen_models['PbergheiANKA_without_ortho'] = cobra.io.read_sbml_model('gf_no_ortho_PbergheiANKA.xml')
essentiality_screen_models['PcynomolgiB_without_ortho'] = cobra.io.read_sbml_model('gf_no_ortho_PcynomolgiB.xml')
essentiality_screen_models['PvivaxSal1_without_ortho'] = cobra.io.read_sbml_model('gf_no_ortho_PvivaxSal1.xml')
essentiality_screen_models['PknowlesiH_without_ortho'] = cobra.io.read_sbml_model('gf_no_ortho_PknowlesiH.xml')
essentiality_screen_models['PcynomolgiB_without_ortho'] = cobra.io.read_sbml_model('gf_no_ortho_PcynomolgiB.xml')

os.chdir("/home/mac9jc/paradigm/data/published_models")
essentiality_screen_models['pfal2018'] = cobra.io.read_sbml_model('pfal2018_abdel_haleem.xml')
essentiality_screen_models['pviv2018'] = cobra.io.read_sbml_model('pviv2018_abdel_haleem.xml')
essentiality_screen_models['pber2018'] = cobra.io.read_sbml_model('pber2018_abdel_haleem.xml')
essentiality_screen_models['pkno2018'] = cobra.io.read_sbml_model('pkno2018_abdel_haleem.xml')
essentiality_screen_models['pcyn2018'] = cobra.io.read_sbml_model('pcyn2018_abdel_haleem.xml')
essentiality_screen_models['ipfa2017'] = cobra.io.read_sbml_model('ipfa2017_chiappino_pepe.xml')
essentiality_screen_models['tg2015'] = cobra.io.read_sbml_model('tg2015_tymoshenko.xml')

os.chdir("/home/mac9jc/paradigm/models")
essentiality_screen_models['iPfal19'] = cobra.io.read_sbml_model('iPfal19.xml')

gene_essentiality_screen_results_raw= dict()
gene_essentiality_screen_results_interpreted = dict()

for species, model in essentiality_screen_models.items():
    raw_results = dict()
    interpreted_results = dict()

    print(species+', gene essenitality screen')
    
    use_second_biomass = False
    if species == 'tg2015':
        model.objective = "Biomass"
        model.reactions.get_by_id('Biomass').upper_bound = 1000.
        model.reactions.get_by_id('Biomass').lower_bound = 0.
    elif species == 'ipfa2017':
        model.objective == 'Biomass_rxn_c'
        model.reactions.get_by_id('Biomass_rxn_c').upper_bound = 1000.
       	model.reactions.get_by_id('Biomass_rxn_c').lower_bound = 0.
    elif species in ['iPfal19','pfal2018','pviv2018','pber2018','pkno2018','pcyn2018']:
        model.objective = "biomass"
        model.reactions.get_by_id('biomass').upper_bound = 1000.
       	model.reactions.get_by_id('biomass').lower_bound = 0.
    elif species in ['TgondiiGT1','TgondiiME49','ChominisTU502_2012','CparvumIowaII']:
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
        pd.DataFrame.from_dict(gene_essentiality_screen_results_raw[species], orient='index').to_csv("/home/mac9jc/paradigm/data/results/gene_essentiality/gene_essentiality_matrix_{}_{}.csv".format(species,day))
        pd.DataFrame.from_dict(gene_essentiality_screen_results_interpreted[species], orient='index').to_csv("/home/mac9jc/paradigm/data/results/gene_essentiality/gene_essentiality_matrix_interpreted_{}_{}.csv".format(species,day))

    else:
       	print('wroking on second biomass reactions')
        gene_essentiality_screen_results_raw[species+'_generic_biomass'] = raw_results
        gene_essentiality_screen_results_interpreted[species+'_generic_biomass'] = interpreted_results
        
        pd.DataFrame.from_dict(gene_essentiality_screen_results_raw[species+'_generic_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/results/gene_essentiality/gene_essentiality_matrix_{}_{}.csv".format(species+'_generic_biomass',day))
        pd.DataFrame.from_dict(gene_essentiality_screen_results_interpreted[species+'_generic_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/results/gene_essentiality/gene_essentiality_matrix_interpreted_{}_{}.csv".format(species+'_generic_biomass',day))

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

        pd.DataFrame.from_dict(gene_essentiality_screen_results_raw[species+'_species_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/results/gene_essentiality/gene_essentiality_matrix_{}_{}.csv".format(species+'_species_biomass',day))
        pd.DataFrame.from_dict(gene_essentiality_screen_results_interpreted[species+'_species_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/results/gene_essentiality/gene_essentiality_matrix_interpreted_{}_{}.csv".format(species+'_species_biomass',day))



