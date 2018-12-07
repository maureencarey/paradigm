
## LOAD MODELS

import cobra
import pandas as pd
import os
import subprocess
import glob
import json
from cobra import Model, Reaction, Metabolite
import helper_functions as hf

# matrix of reaction presence
dict_o_lists = dict()
for model_name, model in model_dict.items():
    dict_o_lists[model_name] = [x.id for x in model.reactions] # get all reaction IDs
unique_elements = set().union(*dict_o_lists.values()) # find all unique reeaction IDs

df = pd.DataFrame(index = unique_elements,columns = dict_o_lists.keys()) # save as matrix/df
for species,list_genes in dict_o_lists.items():
    for x in unique_elements:
        if x in list_genes:
            y = 1
        else:
            y = 0
        df[species][x] = y

os.chdir(data_path)
df.to_csv('reaction_matrix_july25.csv', header=True, index=True)

# plot
proc=subprocess.Popen(["Rscript", "/Users/maureencarey/local_documents/work/comparative_parasite_models/code/upsetR_models_uncurated.R"])
# generating upsetR_all_models_OG_carvemeMAY25.png and others


