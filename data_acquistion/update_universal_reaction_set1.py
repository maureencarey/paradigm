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
logging.basicConfig(filename='update_universal_round1.log', level=logging.INFO, filemode='w')
logger = logging.getLogger(__name__)

os.chdir(model_path)
model = cobra.io.load_json_model('universal_model.json')
cobra.manipulation.modify.escape_ID(model)

# remove Biomass reactions
rxn_list_to_delete = [r.id for r in model.reactions if r.id.startswith('BIOMASS_')]
model.remove_reactions(rxn_list_to_delete)
rxn_list_to_delete = [r.id for r in model.reactions if 'biomass' in r.id]
model.remove_reactions(rxn_list_to_delete)
rxn_list_to_delete = [r.id for r in model.reactions if 'bof' in r.id]
model.remove_reactions(rxn_list_to_delete)
logger.info('removed biomasses from universal')

# add exchange reactions
for met in model.metabolites:
    if met.id.endswith('_e'):
        if 'EX_'+met.id not in model.reactions:
            model.add_boundary(met,type = "exchange")

# add full met and rxn info
met_counter = 0
rxn_counter = 0
for met in model.metabolites:
    if met_counter % 100 == 0:
        logger.info(met_counter)
    met_counter = met_counter + 1
    model = hf3.add_full_met_info(model, met, hf3.met_ids_without_comp(model,met.id))
for rxn in model.reactions:
    if rxn_counter % 100 == 0:
        logger.info(rxn_counter)
    rxn_counter = rxn_counter +1
    model = hf3.add_full_rxn_info(model, rxn, rxn.id)

# fix charges and formulas
model = hf3.fix_charge_or_formula(model)

# Add SBO terms
model = hf3.add_sbo_terms(model)
logger.info('fixed SBO terms')

## remove reactions that have met ids that are specifically human # THIS HAS TO BE DONE PRIOR TO ADDING iPFAL REACTIONS bc we want to keep SM_host and SMPD3l_host
remove_r = list()
add_r = list()
for rxn in model.reactions:
    if len([m.id for m in rxn.metabolites if '_hs_' in m.id])>0:
        remove_r.append(rxn)
        new_rxn = rxn.copy()
        new_rxn.id = rxn.id.replace('_hs_','_')+'nhs'
        new_rxn.name = rxn.name.replace('homo sapiens',',')+' not specific to homo sapiens'
        for met in new_rxn.metabolites:
            if '_hs_' in met.id:
                new_met = met.copy()
                new_met.id = met.id.replace('_hs_', '_')
                new_rxn.metabolites[new_met] = new_rxn.metabolites[met]
                new_rxn.metabolites[met] = 0.
        add_r.append(new_rxn)
add_r = list(set(add_r))
remove_r = list(set(remove_r))
model.add_reactions(add_r)
model.remove_reactions(remove_r)

os.chdir(model_path)
cobra.io.save_json_model(model,  'universal_model_updated.json')
cobra.io.write_sbml_model(model,  'universal_model_updated.xml')

