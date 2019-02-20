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
#data_path = "/Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/data"
#model_path = "/Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/models"

os.chdir(model_path)
universal = cobra.io.load_json_model('universal_model_oct26_2018.json')

cobra.manipulation.modify.escape_ID(universal)
logger.info('finished escape')

def met_ids_without_comp(met_id):
    # only one id listed
    # print list of metabolites without the compartment associated
    # this needs to be updated if you have different compartments than those listed below
    m = met_id
    if m.endswith('_c') or m.endswith('_e') or m.endswith('_f') or \
    m.endswith('_g') or m.endswith('_h') or m.endswith('_i') or \
    m.endswith('_l') or m.endswith('_m') or m.endswith('_n') or \
    m.endswith('_p') or m.endswith('_r') or m.endswith('_s') or \
    m.endswith('_u') or m.endswith('_v') or m.endswith('_x'):
        id_withou_c = m[:-2]
    elif m.endswith('_cx') or m.endswith('_um') or m.endswith('_im') \
    or m.endswith('_ap') or m.endswith('_fv') or m.endswith('_cm'):
        id_withou_c = m[:-3]
    else:
        print('unknown compartment')
        print(m)
        id_withou_c = ''

    return(id_withou_c)

# remove reactions
rxn_list_to_delete = [r.id for r in universal.reactions if r.id.startswith('BIOMASS_')]
universal.remove_reactions(rxn_list_to_delete)
logger.info('removed extra biomass reactions')

# add full met info
for met in universal.metabolites:

    met_id = met_ids_without_comp(met.id)

    m = requests.get('http://bigg.ucsd.edu/api/v2/universal/metabolites/{}'.format(met_id))
    x = m.json()

    if met.name == '':
        met.name = x['name']

    # fix current reaction.annotation field that is in list form
    # this fix will allow universal to be written as as a xml file, not just json
    # however save the info just in case it is not duplicated elsewhere

    # rxn.notes are currently {'original_bigg_id':[id_string]}
    annot_list = met.annotation
    if isinstance(annot_list, list):
        if len(annot_list) >1:
            for sub_list in annot_list:
                if sub_list[0] not in met.notes.keys():
                    met.notes[sub_list[0]] = sub_list[1]
                else:
                    met.notes[sub_list[0]] = [met.notes[sub_list[0]]].append(sub_list[1])
        else:
            met.notes[sub_list[0]] = sub_list[0]

    met.annotation = x['database_links']
    logger.info(met.id)
    logger.info(x['database_links'])
    if x['formulae'] != [] and len(x['formulae'])>0:
        met.formula = x['formulae'][0]
    if x['charges'] != [] and len(x['charges'])>0:
        met.charge = x['charges'][0]

# add full reaction info
need_info = list()
for rxn in universal.reactions:

    rxn_id = rxn.id
    m = requests.get('http://bigg.ucsd.edu/api/v2/universal/reactions/{}'.format(rxn_id))
    x = m.json()

    if rxn.name == '':
        rxn.name = x['name']

    # fix current reaction.annotation field that is in list form
    # this fix will allow universal to be written as as a xml file, not just json
    # however save the info just in case it is not duplicated elsewhere

    # rxn.notes are currently {'original_bigg_id':[id_string]}
    annot_list = rxn.annotation
    if isinstance(annot_list, list):
        if len(annot_list) >1:
            for sub_list in annot_list:
                if sub_list[0] not in rxn.notes.keys():
                    rxn.notes[sub_list[0]] = sub_list[1]
                else:
                    rxn.notes[sub_list[0]] = [rxn.notes[sub_list[0]]].append(sub_list[1])
        else:
            rxn.notes[sub_list[0]] = sub_list[0]
    if rxn.reaction == '':
        rxn.reaction = x['reaction_string']
    rxn.annotation = x['database_links']
    logger.info(rxn.id)
    logger.info(x['database_links'])
    # also have info on x['pseudoreaction']
    # COULD ADD SBO TERM HERE

for met in universal.metabolites:
    if isinstance(met.charge, list):
        if len(met.charge) == 0:
            met.charge = int(0)
        else:
            met.charge = met.charge[0]
    if isinstance(met.formula, list):
        if len(met.formula) == 0:
            met.formula = ''
        else:
            met.formula = met.formula[0]


universal.reactions.PGM.lower_bound = -1000. # make reversible

# SAVE
os.chdir(model_path)
cobra.io.save_json_model(universal, 'universal_model_updated.json')
cobra.io.write_sbml_model(universal, 'universal_model_updated.xml')

#universal = cobra.io.load_json_model('universal_model_updated.json')
rxn_list = [r.id for r in universal.reactions]
met_list = [m.id for m in universal.metabolites]

# same for Pf model
os.chdir(model_path)
model = cobra.io.load_json_model('iPfal18.json')
cobra.manipulation.modify.escape_ID(model)

# add full met info
for met in model.metabolites:
    met_id = met_ids_without_comp(met.id)

    if met_id+'_c' in met_list:
        
        m = requests.get('http://bigg.ucsd.edu/api/v2/universal/metabolites/{}'.format(met_id))
        x = m.json()
        if met.name == '':
            met.name = x['name']
        
        # fix current reaction.annotation field that is in list form
        # this fix will allow universal to be written as as a xml file, not just json
        # however save the info just in case it is not duplicated elsewhere
        
        # rxn.notes are currently {'original_bigg_id':[id_string]}
        annot_list = met.annotation
        if isinstance(annot_list, list):
            if len(annot_list) >1:
                for sub_list in annot_list:
                    if sub_list[0] not in met.notes.keys():
                        met.notes[sub_list[0]] = sub_list[1]
                    else:
                        met.notes[sub_list[0]] = [met.notes[sub_list[0]]].append(sub_list[1])
            else:
                met.notes[sub_list[0]] = sub_list[0]

        met.annotation = x['database_links'] # check if compartment info

        logger.info(met.id)
        logger.info(x['database_links'])
        if x['formulae'] != [] and len(x['formulae'])>0:
            met.formula = x['formulae'][0]
        if x['charges'] != [] and len(x['charges'])>0:
            met.charge = x['charges'][0]

# add full reaction info
need_info = list()
for rxn in model.reactions:
    
    rxn_id = rxn.id

    if rxn_id in rxn_list:
        m = requests.get('http://bigg.ucsd.edu/api/v2/universal/reactions/{}'.format(rxn_id))
        x = m.json()
        if rxn.name == '':
            rxn.name = x['name']
        
        # fix current reaction.annotation field that is in list form
        # this fix will allow universal to be written as as a xml file, not just json
        # however save the info just in case it is not duplicated elsewhere
        
        # rxn.notes are currently {'original_bigg_id':[id_string]}
        annot_list = rxn.annotation
        if isinstance(annot_list, list):
            if len(annot_list) >1:
                for sub_list in annot_list:
                    if sub_list[0] not in rxn.notes.keys():
                        rxn.notes[sub_list[0]] = sub_list[1]
                    else:
                        rxn.notes[sub_list[0]] = [rxn.notes[sub_list[0]]].append(sub_list[1])
            else:
                rxn.notes[sub_list[0]] = sub_list[0]
        rxn.annotation = x['database_links']
        logger.info(rxn.id)
        logger.info(x['database_links'])
        # also have info on x['pseudoreaction']
        # COULD ADD SBO TERM HERE
# also have info on x['pseudoreaction']

model.reactions.get_by_id('DM_biomass_c').name = 'unblock biomass'
model.reactions.get_by_id('DM_biomass_c').id = 'DM_bm'

for met in model.metabolites:
    if isinstance(met.charge, list):
        if len(met.charge) == 0:
            met.charge = int(0)
        else:
            met.charge = met.charge[0]
    if isinstance(met.formula, list):
        if len(met.formula) == 0:
            met.formula = ''
        else:
            met.formula = met.formula[0]
    if met.compartment == 'food vacuole':
        met.compartment = 'food_vacuole'


## add notes field from previous curation rounds

os.chdir(data_path)
curation_file1 = pd.read_csv("")
curation_file2 = pd.read_csv("")

for rxn in model.reactions:
    if rxn.id in curation_file1['reaction']:
    if rxn.id in curation_file2['reaction']:
    
    if isinstance(rxn.notes,dict):
        rxn.notes['notes'] =
        rxn.notes['reference'] =
    elif isinstance(rxn.notes,list):
    else:



cobra.io.save_json_model(model,  'iPfal18_updated.json')
cobra.io.write_sbml_model(model,  'iPfal18_updated.xml')

