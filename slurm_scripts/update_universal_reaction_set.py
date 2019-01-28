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

#data_path = "/home/mac9jc/paradigm/data/"
#model_path = "/home/mac9jc/paradigm/models/"
data_path = "/Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/data"
model_path = "/Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/models"

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
    m = met_id
    if m.endswith('_c') or m.endswith('_e') or m.endswith('_f') or \
    m.endswith('_g') or m.endswith('_h') or m.endswith('_i') or \
    m.endswith('_l') or m.endswith('_m') or m.endswith('_n') or \
    m.endswith('_p') or m.endswith('_r') or m.endswith('_s') or \
    m.endswith('_u') or m.endswith('_v') or m.endswith('_x'):
        id_withou_c = m[:-2]
    elif m.endswith('_cx') or m.endswith('_um') or m.endswith('_im') \
    or m.endswith('_ap') or m.endswith('_fv'):
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
    
    # fix current reaction.annotations field that is in list form
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

    met.annotations = x['database_links']
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

    # fix current reaction.annotations field that is in list form
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
    rxn.reaction = x['reaction_string']
    rxn.annotation = x['database_links']
    # also have info on x['pseudoreaction']
    # COULD ADD SBO TERM HERE

# SAVE
os.chdir(model_path)
cobra.io.save_json_model(universal, 'universal_model_updated.json')
cobra.io.write_sbml_model(universal, 'universal_model_updated.xml')

#rxn_list = [r.id for r in universal.reactions]
#met_list = [m.id for m in universal.metabolites]
#
## same for Pf model
#os.chdir(model_path)
#model = cobra.io.load_json_model('iPfal18.json')
#
## add full met info
#for met in model.metabolites:
#    met_id = met_ids_without_comp(met.id)
#
#    if met_id+'_c' in met_list:
#        m = requests.get('http://bigg.ucsd.edu/api/v2/universal/metabolites/{}'.format(met_id+'_c'))
#        x = m.json()
#        if met.name == '':
#            met.name = x['name']
#        if met.notes = '':
#            met.notes = met.annotation
#        met.annotation = x['database_links'] # check if compartment info
        # if x['formulae'] is not []:
    #        met.formula = x['formulae'][0]
#        if x['charges'] is not []:
    #        met.charge = x['charges'][0]
#
## add full reaction info
#need_info = list()
#for rxn in model.reactions:
#    rxn_id = rxn.id
#
#    if rxn_id in rxn_list:
#        m = requests.get('http://bigg.ucsd.edu/api/v2/universal/reactions/{}'.format(rxn_id))
#        x = m.json()
#        if rxn.name == '':
#            rxn.name = x['name']
#        if rxn.notes = '':
#            rxn.notes = rxn.annotation
#        rxn.annotation = x['database_links'] # check if compartment info
## also have info on x['pseudoreaction']
#
# model.reactions.get_by_id('biomass_s').name = 'biomass sink'
# model.reactions.get_by_id('biomass_s').id = 'SK_bm'
#save_json_model(model,  'iPfal18_updated.json')

