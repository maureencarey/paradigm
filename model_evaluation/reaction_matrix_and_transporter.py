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
presence_matrix_of_transporters.to_csv("/home/mac9jc/paradigm/data/results/transporter_presence_before_gapfilling_jan.csv")

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
presence_matrix_of_reactions.to_csv("/home/mac9jc/paradigm/data/results/rxn_presence_before_gapfilling_jan.csv")
