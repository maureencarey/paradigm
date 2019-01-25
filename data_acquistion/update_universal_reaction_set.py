import cobra
import os
import pandas as pd
import numpy as np
import glob
from cobra.core import Gene, Metabolite, Reaction
from datetime import date
from datetime import datetime
import requests
import logging

day = datetime.now().strftime('%d_%m_%Y')
logging.basicConfig(filename='update_reaction_set_{}.log'.format(day), level=logging.INFO, filemode='w')
logger = logging.getLogger(__name__)

logger.info('begin: update of universal reaction set')

data_path = "/home/mac9jc/paradigm/data/"
model_path = "/home/mac9jc/paradigm/models/"

os.chdir(model_path)
universal = cobra.io.load_json_model('universal_model_oct26_2018.json')

os.chdir(data_path)
met_document = pd.read_table('bigg_metabolites.txt')
logger.info('loaded universal model and bigg met info')

cobra.manipulation.modify.escape_ID(universal)
logger.info('finished escape')

def met_ids_without_comp(met_id):
    # only one id listed
    # print list of metabolites without the compartment associated
    # this needs to be updated if you have different compartments than those listed below

    if m.id.endswith('_c') or m.id.endswith('_e') or m.id.endswith('_f') or \
    m.id.endswith('_g') or m.id.endswith('_h') or m.id.endswith('_i') or \
    m.id.endswith('_l') or m.id.endswith('_m') or m.id.endswith('_n') or \
    m.id.endswith('_p') or m.id.endswith('_r') or m.id.endswith('_s') or \
    m.id.endswith('_u') or m.id.endswith('_v') or m.id.endswith('_x'):
        id_withou_c = m.id[:-2]
    elif m.id.endswith('_cx') or m.id.endswith('_um') or m.id.endswith('_im') \
    or m.id.endswith('_ap') or m.id.endswith('_fv'):
        id_withou_c = m.id[:-3]
    else:
        print('unknown compartment')
        print(m.id)
        id_withou_c = ''

    return(id_withou_c)

# remove reactions
rxn_list_to_delete = [r.id for r in universal.reactions if r.id.startswith('BIOMASS_')]
universal.remove_reactions(rxn_list_to_delete)
logger.info('removed extra biomass reactions')

# add full met info
for met in universal.metabolites:
    met_id = met_ids_without_comp(met.id)
    m = requests.get('http://bigg.ucsd.edu/api/v2/universal/metabolites/{}'.format(met_id)))
    x = m.json()
    met.name = x['name']
    met.notes = x['database_links']
    met.formula = x['formulae'][0]
    met.charge = x['charges'][0]

# add full reaction info
need_info = list()
for rxn in universal.reactions:
    rxn_id = rxn.id
    m = requests.get('http://bigg.ucsd.edu/api/v2/universal/reactions/{}'.format(rxn_id))
    x = m.json()
    rxn.name = x['name']
    rxn.reaction = x['reaction_string']
    rxn.name = x['name']
    rxn.notes = x['database_links']
    # also have info on x['pseudoreaction']

### WHAT ABOUT METS AND REACTIONS IN OTHER COMPARTMENTS

## SAVE MODEL


