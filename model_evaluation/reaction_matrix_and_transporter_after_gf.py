import cobra
import os
import pandas as pd
import numpy as np
import glob
from cobra.core import Gene, Metabolite, Reaction
from datetime import datetime

day = datetime.now().strftime('%d_%m_%Y')

model_dict = dict()
path = "/home/mac9jc/paradigm/models/"
os.chdir(path)
for filename in glob.glob(os.path.join(path, 'gf_*.xml')):
    key = filename.split('/')[len(filename.split('/'))-1]
    key = key[:-4]
    key = key[3:]
    print(key)
    model_dict[key] = cobra.io.read_sbml_model(filename)

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
presence_matrix_of_reactions.to_csv("/home/mac9jc/paradigm/data/results/rxn_presence_after_gapfilling_{}.csv".format(day))
