import cobra
import os
import pandas as pd
import numpy as np
import glob
from cobra.core import Gene, Metabolite, Reaction

## ##### get essential reactions
os.chdir("/home/mac9jc/paradigm/models")
essentiality_screen_models = dict()
essentiality_screen_models['TgondiiGT1'] = cobra.io.load_json_model('gf_TgondiiGT1.json')
essentiality_screen_models['TgondiiME49'] = cobra.io.load_json_model('gf_TgondiiME49.json')
essentiality_screen_models['Pfalciparum3D7'] = cobra.io.load_json_model('gf_Pfalciparum3D7.json')
essentiality_screen_models['PbergheiANKA'] = cobra.io.load_json_model('gf_PbergheiANKA.json')
essentiality_screen_models['PcynomolgiB'] = cobra.io.load_json_model('gf_PcynomolgiB.json')
essentiality_screen_models['PvivaxSal1'] = cobra.io.load_json_model('gf_PvivaxSal1.json')
essentiality_screen_models['ChominisTU502_2012'] = cobra.io.load_json_model('gf_ChominisTU502_2012.json')
essentiality_screen_models['CparvumIowaII'] = cobra.io.load_json_model('gf_CparvumIowaII.json')
essentiality_screen_models['PknowlesiH'] = cobra.io.load_json_model('gf_PknowlesiH.json')
essentiality_screen_models['PcynomolgiB'] = cobra.io.load_json_model('gf_PcynomolgiB.json')

print('model set 1 loaded')

essentiality_screen_models['Pfalciparum3D7_without_ortho'] = cobra.io.load_json_model('gf_without_ortho_Pfalciparum3D7.json')
essentiality_screen_models['PbergheiANKA_without_ortho'] = cobra.io.load_json_model('gf_without_ortho_PbergheiANKA.json')
essentiality_screen_models['PcynomolgiB_without_ortho'] = cobra.io.load_json_model('gf_without_ortho_PcynomolgiB.json')
essentiality_screen_models['PvivaxSal1_without_ortho'] = cobra.io.load_json_model('gf_without_ortho_PvivaxSal1.json')
essentiality_screen_models['PknowlesiH_without_ortho'] = cobra.io.load_json_model('gf_without_ortho_PknowlesiH.json')
essentiality_screen_models['PcynomolgiB_without_ortho'] = cobra.io.load_json_model('gf_without_ortho_PcynomolgiB.json')

print('model set 2 loaded')

os.chdir("/home/mac9jc/paradigm/data/published_models")
essentiality_screen_models['pfal2018'] = cobra.io.read_sbml_model('pfal2018_abdel_haleem.xml')
essentiality_screen_models['pviv2018'] = cobra.io.read_sbml_model('pviv2018_abdel_haleem.xml')
essentiality_screen_models['pber2018'] = cobra.io.read_sbml_model('pber2018_abdel_haleem.xml')
essentiality_screen_models['pkno2018'] = cobra.io.read_sbml_model('pkno2018_abdel_haleem.xml')
essentiality_screen_models['pcyn2018'] = cobra.io.read_sbml_model('pcyn2018_abdel_haleem.xml')
essentiality_screen_models['ipfa2017'] = cobra.io.read_sbml_model('ipfa2017_chiappino_pepe.xml')
essentiality_screen_models['tg2015'] = cobra.io.read_sbml_model('tg2015_tymoshenko.xml')

print('model set 3 loaded')

os.chdir("/home/mac9jc/paradigm/models")
essentiality_screen_models['iPfal18'] = cobra.io.load_json_model('iPfal18_updated.json')

essentiality_screen_results_raw= dict()
essentiality_screen_results_interpreted = dict()

print('model set 4 loaded')

for species, model in essentiality_screen_models.items():
    raw_results = dict()
    interpreted_results = dict()
    
    use_second_biomass = False
    
    print(species+', rxn essenitality screen')
    
    use_second_biomass = False
    if species == 'tg2015':
        model.objective = "Biomass"
        model.reactions.get_by_id('Biomass').upper_bound = 1000.
       	model.reactions.get_by_id('Biomass').lower_bound = 0.
    elif species == 'ipfa2017':
        model.objective == 'Biomass_rxn_c'
        model.reactions.get_by_id('Biomass_rxn_c').upper_bound = 1000.
       	model.reactions.get_by_id('Biomass_rxn_c').lower_bound = 0.
    elif species in ['iPfal18','pfal2018','pviv2018','pber2018','pkno2018','pcyn2018']:
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

    print('set biomass')

    max_biomass = model.slim_optimize()
    if max_biomass < 0.01:
        print('model doesnt grow')
        print(species)

    # knockout and record growth
    for rxn in model.reactions:
        if rxn.id == 'generic_biomass' or rxn.id == 'biomass':
            continue
        with model as cobra_model:
            cobra_model.reactions.get_by_id(rxn.id).knock_out()
            f = cobra_model.slim_optimize()
            if max_biomass < 0.01:
                print(species)
            else:
                if f < 0.1*max_biomass:
                    interpreted_results[rxn.id] = 'lethal'
                else:
                    interpreted_results[rxn.id] = 'nonlethal'
                raw_results[rxn.id] = f/max_biomass

    # save
    if not use_second_biomass:
        
        essentiality_screen_results_raw[species] = raw_results
        essentiality_screen_results_interpreted[species] = interpreted_results
        cols = ['gene_id', 'normalized_growth']
        pd.DataFrame.from_dict(essentiality_screen_results_raw[species], orient='index').to_csv("/home/mac9jc/paradigm/data/rxn_essentiality_matrix_{}.csv".format(species))
        pd.DataFrame.from_dict(essentiality_screen_results_interpreted[species], orient='index').to_csv("/home/mac9jc/paradigm/data/rxn_essentiality_matrix_interpreted_{}.csv".format(species))
        
    else:
        essentiality_screen_results_raw[species+'_generic_biomass'] = raw_results
        essentiality_screen_results_interpreted[species+'_generic_biomass'] = interpreted_results
        cols = ['gene_id', 'normalized_growth']
        pd.DataFrame.from_dict(essentiality_screen_results_raw[species+'_generic_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/rxn_essentiality_matrix_{}.csv".format(species+'_generic_biomass'))
        pd.DataFrame.from_dict(essentiality_screen_results_interpreted[species+'_generic_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/rxn_essentiality_matrix_interpreted_{}.csv".format(species+'_generic_biomass'))

        
        print(species+', PLASMO SPECIFIC rxn essenitality screen')
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
        for rxn in model.reactions:
            if rxn.id == 'generic_biomass' or rxn.id == 'biomass':
                continue
            with model as cobra_model:
                cobra_model.reactions.get_by_id(rxn.id).knock_out()
                f = cobra_model.slim_optimize()
                if max_biomass < 0.01:
       	            print(species)
                else:
                    if f < 0.1*max_biomass:
                        interpreted_results[rxn.id] = 'lethal'
                    else:
                        interpreted_results[rxn.id] = 'nonlethal'
                    raw_results[rxn.id] = f/max_biomass

        # save
        essentiality_screen_results_raw[species+'_species_biomass'] = raw_results
        essentiality_screen_results_interpreted[species+'_species_biomass'] = interpreted_results

        pd.DataFrame.from_dict(essentiality_screen_results_raw[species+'_species_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/rxn_essentiality_matrix_{}.csv".format(species+'_species_biomass'))
        pd.DataFrame.from_dict(essentiality_screen_results_interpreted[species+'_species_biomass'], orient='index').to_csv("/home/mac9jc/paradigm/data/rxn_essentiality_matrix_interpreted_{}.csv".format(species+'_species_biomass'))

list_o_reactions2 = list()
for species, model in essentiality_screen_models.items():
    list_o_reactions2.append([x.id for x in model.reactions])
list_o_reactions2 = [val for sublist in list_o_reactions2 for val in sublist]
list_o_reactions2 = list(set(list_o_reactions2))

matrix_of_essentiality = pd.DataFrame(index = list_o_reactions2,columns=essentiality_screen_results_raw.keys())
matrix_interpreted = pd.DataFrame(index = list_o_reactions2,columns=essentiality_screen_results_interpreted.keys())
for species_long, rxn_list in essentiality_screen_results_raw.items():
    print(species_long)
    if '_generic' in species_long:
        species = species_long.split('_generic')[0]
    elif '_species' in species_long:
        species = species_long.split('_species')[0]
    else:
        print('ERROR, what biomass are we using?')
    for rxn in list_o_reactions2:
        if rxn in essentiality_screen_results_raw[species_long].keys():
            matrix_of_essentiality.loc[rxn,species_long] = essentiality_screen_results_raw[species_long][rxn]
            matrix_interpreted.loc[rxn,species_long] = essentiality_screen_results_interpreted[species_long][rxn]
        else:
            matrix_of_essentiality.loc[rxn,species_long] = 'NA'
            matrix_interpreted.loc[rxn,species_long] = 'NA'
matrix_of_essentiality.to_csv("/home/mac9jc/paradigm/data/rxn_essentiality_matrix_jan.csv")
matrix_interpreted.to_csv("/home/mac9jc/paradigm/data/rxn_essentiality_matrix_interpreted_jan.csv")



