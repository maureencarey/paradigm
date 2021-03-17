import cobra
import os
import pandas as pd
from cobra.core import Gene, Metabolite, Reaction
import requests
import time
import logging
import sys
sys.path.append(os.path.abspath("/home/mac9jc/paradigm/"))
import helper_functions as hf

model_path = "/home/mac9jc/paradigm/models/"
log_path = "/home/mac9jc/paradigm/model_generation_logs/"

os.chdir(log_path)
logging.basicConfig(filename='update_universal_round1.log', level=logging.INFO, filemode='w')
logger = logging.getLogger(__name__)

os.chdir(model_path)
model = cobra.io.load_json_model('universal/universal_model.json')
cobra.manipulation.modify.escape_ID(model)
model.repair()

# remove Biomass reactions
rxn_list_to_delete = [r.id for r in model.reactions if r.id.startswith('BIOMASS_')]
model.remove_reactions(rxn_list_to_delete)
rxn_list_to_delete = [r.id for r in model.reactions if 'biomass' in r.id]
model.remove_reactions(rxn_list_to_delete)
rxn_list_to_delete = [r.id for r in model.reactions if 'bof' in r.id]
model.remove_reactions(rxn_list_to_delete)
logger.info('removed biomasses from universal')
model.repair()

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
    model = hf.add_full_met_info(model, met, hf.met_ids_without_comp(model,met.id))
for rxn in model.reactions:
    if rxn_counter % 100 == 0:
        logger.info(rxn_counter)
    rxn_counter = rxn_counter +1
    model = hf.add_full_rxn_info(model, rxn, rxn.id)

# fix charges and formulas
model = hf.fix_charge_or_formula(model)

# Add SBO terms
model = hf.add_sbo_terms(model)
logger.info('fixed SBO terms')
model.repair()

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
model.repair()

compartment_dict = {'c': 'cytoplasm', 'e': 'extracellular','m': 'mitochondrion', 'fv': 'food vacuole', 'ap':'apicoplast','k':'kinetoplast','glc':'glycosome', 'p':'periplasm','x':'peroxisome/glyoxysome','r':'endoplasmic reticulum','v':'vacuole','n':'nucleus','g':'golgi apparatus','u':'thylakoid','l':'lysosome','h':'chloroplast','f':'flagellum','s':'eyespot','im':'intermembrane space of mitochondria','cx':'carboxyzome','um':'thylakoid membrane','cm':'cytosolic membrane','i':'not_provided_by_bigg'}

for met in model.metabolites:
    comp = hf.get_comp(model,met.id)
    comp = comp[1:]
    if comp in compartment_dict:
        comp_use = compartment_dict[comp]
    else: comp_use = 'no compartment'
    model.metabolites.get_by_id(met.id).compartment = comp_use
model.repair()

os.chdir(model_path)
cobra.io.save_json_model(model,  './universal/universal_model_updated.json')
cobra.io.write_sbml_model(model,  './universal/universal_model_updated.xml')

