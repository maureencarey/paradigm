import cobra
import os
import pandas as pd
from cobra.core import Gene, Metabolite, Reaction
import requests
import time
import logging
import helper_functions_3 as hf3

model_path = "/home/mac9jc/paradigm/models/"
log_path = "/home/mac9jc/paradigm/model_generation_logs/"

os.chdir(log_path)
logging.basicConfig(filename='update_universal_round2.log', level=logging.INFO, filemode='w')
logger = logging.getLogger(__name__)

os.chdir(model_path)
model = cobra.io.load_json_model('universal_model_updated_r1.json')
cobra.manipulation.modify.escape_ID(model)

os.chdir(model_path)
model_dict = dict()
for filename in glob.glob(os.path.join(path, 'final_denovo_*.json')):
    key = filename.split('/')[len(filename.split('/'))-1]
    key = key[:-5]
    key = key[13:]
    logger.info(key)
    model_dict[key] = cobra.io.load_json_model(filename)

for species, denovo_mod in model_dict.items():
    for reaction in denovo_mod.reactions:
        if reaction.id not in [r.id for r in model.reactions]:
            model.add_reactions([reaction.copy()])

# remove Biomass reactions
rxn_list_to_delete = [r.id for r in model.reactions if r.id.startswith('BIOMASS_')]
model.remove_reactions(rxn_list_to_delete)
rxn_list_to_delete = [r.id for r in model.reactions if 'biomass' in r.id]
model.remove_reactions(rxn_list_to_delete)
logger.info('removed biomasses from universal')

# add exchange reactions
for met in model.metabolites:
    if met.id.endswith('_e'):
        if 'EX_'+met.id not in model.reactions:
            model.add_boundary(met,type = "exchange")

# Add SBO terms
model = hf3.add_sbo_terms(model)
logger.info('fixed SBO terms')

os.chdir(model_path)
cobra.io.save_json_model(model,  'universal_model_updated_for_gapfilling.json')
cobra.io.write_sbml_model(model,  'universal_model_updated_for_gapfilling.xml')

