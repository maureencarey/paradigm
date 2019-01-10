import cobra
import os
import pandas as pd
import numpy as np
import glob
from cobra.core import Gene, Metabolite, Reaction

model_dict = dict()
path = "/home/mac9jc/paradigm/models/"
os.chdir(path)
for filename in glob.glob(os.path.join(path, 'final_denovo_*.json')):
    key = filename.split('/')[len(filename.split('/'))-1]
    key = key[:-5]
    key = key[13:]
#    print(key)
    model_dict[key] = cobra.io.load_json_model(filename)
# for filename in glob.glob(os.path.join(path, 'ortho_*.json')):
  #  key = filename.split('/')[len(filename.split('/'))-1]
   # key = key[:-5]
  #  key = key[6:]
  #  print(key)
  #  model_dict[key] = cobra.io.load_json_model(filename)

# these metabolites are transported INTO the cell
imported_mets_dict = dict()
for species, model in model_dict.items():
    reactants = list()
    for rxn in model.reactions:
        reactants.append([x.id for x in rxn.reactants])
    reactants = [val for sublist in reactants for val in sublist]
    transported_mets = list()
    for x in model.metabolites:
        if x.id.endswith('_e') and x.id in reactants:
            transported_mets.append(x.id[:-2])
    imported_mets_dict[species] = list(set(transported_mets))

# all transported mets
met_list = list()
for species, imported_mets in imported_mets_dict.items():
    met_list.append(imported_mets)
met_list = list(set([val for sublist in met_list for val in sublist]))

# get matrix of imported mets
presence_matrix_of_transporters = pd.DataFrame(index = met_list,columns=model_dict.keys())
for species, imported_mets in imported_mets_dict.items():
    for met in met_list:
        if met in imported_mets:
            presence_matrix_of_transporters.loc[met,species] = 1
        else:
            presence_matrix_of_transporters.loc[met,species] = 0
presence_matrix_of_transporters.to_csv("/home/mac9jc/paradigm/data/transporter_presence_before_gapfilling_jan.csv")

# get matrix of reactions
reactions_in_model = dict()
list_o_reactions = list()
for species, model in model_dict.items():
    reactions_in_model[species] = [x.id for x in model.reactions]
    list_o_reactions.append([x.id for x in model.reactions])
list_o_reactions = [val for sublist in list_o_reactions for val in sublist]
list_o_reactions = list(set(list_o_reactions))
presence_matrix_of_reactions = pd.DataFrame(index = list_o_reactions,columns=model_dict.keys())
for species, rxn_list in reactions_in_model.items():
    for rxn in rxn_list:
        if rxn in list_o_reactions:
            presence_matrix_of_reactions.loc[rxn,species] = 1
        else:
            presence_matrix_of_reactions.loc[rxn,species] = 0
presence_matrix_of_reactions.to_csv("/home/mac9jc/paradigm/data/rxn_presence_before_gapfilling_jan.csv")

# get essential reactions
os.chdir("/home/mac9jc/paradigm/models")
essentiality_screen_models = dict()
essentiality_screen_models['TgondiiRH'] = cobra.io.load_json_model('gf_TgondiiRH.json')
essentiality_screen_models['Pfalciparum3D7'] = cobra.io.load_json_model('gf_Pfalciparum3D7.json')
essentiality_screen_models['PbergheiANKA'] = cobra.io.load_json_model('gf_PbergheiANKA.json')
essentiality_screen_models['ChominisTU502_2012'] = cobra.io.load_json_model('gf_ChominisTU502_2012.json')
essentiality_screen_models['CparvumIowaII'] = cobra.io.load_json_model('gf_CparvumIowaII.json')

essentiality_screen_results_raw= dict()
essentiality_screen_results_interpreted = dict()

for species, model in essentiality_screen_models.items():
    raw_results = dict()
    interpreted_results = dict()
    
    # set objective
    model.objective = "generic_biomass"
    
    # don't accidentally use other biomass reaction
    if 'biomass' in [rxn.id for rxn in model.reactions]:
        model.reactions.get_by_id('generic_biomass').upper_bound = 1000.
        model.reactions.get_by_id('generic_biomass').lower_bound = 0.
        model.reactions.get_by_id('biomass').upper_bound = 0.
        model.reactions.get_by_id('biomass').lower_bound = 0.
    
    max_biomass = model.slim_optimize()

    # knockout and record growth
    for rxn in model.reactions:
        if rxn.id == 'generic_biomass' or rxn.id == 'biomass':
            continue
        with model as cobra_model:
            cobra_model.reactions.get_by_id(rxn.id).knock_out()
            f = cobra_model.slim_optimize()
            if f < 0.1*max_biomass:
                interpreted_results[rxn.id] = 'lethal'
            else:
                interpreted_results[rxn.id] = 'nonlethal'
            raw_results[rxn.id] = f

    # save
    essentiality_screen_results_raw[species+'_generic_biomass'] = raw_results
    essentiality_screen_results_interpreted[species+'_generic_biomass'] = interpreted_results

    if species.startswith('P'):
        # set objective
        model.objective = "biomass"
            
        # don't accidentally use other biomass reaction
        model.reactions.get_by_id('biomass').upper_bound = 1000.
        model.reactions.get_by_id('biomass').lower_bound = 0.
        model.reactions.get_by_id('generic_biomass').upper_bound = 0.
        model.reactions.get_by_id('generic_biomass').lower_bound = 0.
            
        max_biomass = model.slim_optimize()
            
        # knockout and record growth
        for rxn in model.reactions:
            if rxn.id == 'generic_biomass' or rxn.id == 'biomass':
                continue
            with model as cobra_model:
                cobra_model.reactions.get_by_id(rxn.id).knock_out()
                f = cobra_model.slim_optimize()
                if f < 0.1*max_biomass:
                    interpreted_results[rxn.id] = 'lethal'
                else:
                    interpreted_results[rxn.id] = 'nonlethal'
                raw_results[rxn.id] = f

        # save
        essentiality_screen_results_raw[species+'_species_biomass'] = raw_results
        essentiality_screen_results_interpreted[species+'_species_biomass'] = interpreted_results

matrix_of_essentiality = pd.DataFrame(index = list_o_reactions,columns=essentiality_screen_results_raw.keys())
for species_long, rxn_list in essentiality_screen_results_raw.items():
    if '_generic' in species:
        species = species_long.split('_generic')[0]
    elif '_species' in species_long:
        species = species_long.split('_species')[0]
    else:
        print('ERROR, what biomass are we using?')
    for rxn in list_o_reactions:
        matrix_of_essentiality.loc[rxn,species] = essentiality_screen_results_raw[species][rxn]
matrix_of_essentiality.to_csv("/home/mac9jc/paradigm/data/rxn_essentiality_matrix_jan.csv")




