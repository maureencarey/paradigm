# INPUT FILENAME???
import cobra
import os
from cobra import Model, Reaction, Metabolite
import pandas as pd
import requests
import argparse


og_path = "/Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm"
data_path = "/Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/data"
model_path = "/Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/models"

os.chdir(data_path)
met_document = pd.read_table('bigg_metabolites.txt')

parser = argparse.ArgumentParser(description='Read in the species model')
parser.add_argument('model_file')
args = parser.parse_args()
model_fname = vars(args)['model_file']

# parse arguments for global variables 
SPECIES_ID = model_fname.split('/')[-1] # ID is model filename minus directory
SPECIES_ID = SPECIES_ID.split('.')[0] # get rid of extension
SPECIES_ID = SPECIES_ID.split('with_biomass_denovo_')[1]

print(SPECIES_ID)
# modified for Rivanna: read in the models
model = cobra.io.load_json_model(model_fname)
print('loaded model')

cobra.manipulation.modify.escape_ID(model)

print('finished escape')

need_info = list()
for met in model.metabolites:
    if met.id in [m.id for m in universal_model.metabolites]:
        if met.id.startswith('protein') or met.id.startswith('lipid'):
            test_string  = 's'
        else:
            met_row = met_document[met_document['bigg_id'] == met.id]
            met_row['universal_bigg_id']
            m = requests.get('http://bigg.ucsd.edu/api/v2/universal/metabolites/{}'.\
            format(met_row['universal_bigg_id'][met_row['universal_bigg_id'].index[0]]))
            x = m.json()
            met.name = x['name']
            met.notes = x['database_links']
            met.formula = x['formulae'][0]
            met.charge = x['charges'][0]
    elif met.id.endswith('_ap') or met.id.endswith('_fv'):
        if met.id[:-3] in [met.id[:-2] for met in universal_model.metabolites]:
            m = requests.get('http://bigg.ucsd.edu/api/v2/universal/metabolites/{}'.\
            format(met.id[:-3]))
            x = m.json()
            met.name = x['name']
            met.notes = x['database_links']
            met.formula = x['formulae'][0]
            met.charge = x['charges'][0]
        else:
            need_info.append(met.id)
    elif met.id[:-2] in [met.id[:-2] for met in universal_model.metabolites]:
        m = requests.get('http://bigg.ucsd.edu/api/v2/universal/metabolites/{}'.\
        format(met.id[:-2]))
        x = m.json()
        met.name = x['name']
        met.notes = x['database_links']
        met.formula = x['formulae'][0]
        met.charge = x['charges'][0]
    else:
        need_info.append(met.id)
print(list(set(need_info)))

## SAVE MODEL
        

       

# FIND/ REMOVE DUPLICATE REACTIONS
# 
# duplicates = dict()
# temp_dict = dict() # get products and reactants for every reaction
# for rxn in pf_model.reactions:
#     rxn_dict = dict()
#     check_rxn_products = rxn.products
#     check_rxn_reactants = rxn.reactants
#     rxn_dict['reactants'] = [x.id for x in check_rxn_reactants]
#     rxn_dict['products'] = [x.id for x in check_rxn_products]
#     temp_dict[rxn.id] = rxn_dict
# 
# print('temp dict made')
# print([key for key in temp_dict.keys() if key not in [x.id for x in pf_model.reactions]])
# 
# for rxn in universal_model.reactions:
#     for key in temp_dict.keys():
#         if key != rxn.id:
#             if [x.id for x in rxn.reactants] == temp_dict[key]['reactants'] and \
#             [x.id for x in rxn.products] == temp_dict[key]['products']:
#                 if rxn.id not in duplicates.keys():
#                     duplicates[rxn.id] = key
#                 elif duplicates[rxn.id] == key or key in duplicates[rxn.id]:
#                     continue
#                 else:
#                      duplicates[rxn.id] = duplicates[rxn.id]+', '+key
#                         
# print('duplicates made')
