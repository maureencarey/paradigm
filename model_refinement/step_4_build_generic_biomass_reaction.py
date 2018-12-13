## Input = final_denovo_SPECIES.json

import cobra
import pandas as pd
import os
import subprocess
import glob
import json
from cobra import Model, Reaction, Metabolite
import helper_functions_2 as hf
import argparse
from datetime import date
from datetime import datetime
import logging

data_path = "/home/mac9jc/paradigm/data"
model_path = "/home/mac9jc/paradigm/models"
os.chdir(model_path)

parser = argparse.ArgumentParser(description='Read in the species model')
parser.add_argument('model_file')
args = parser.parse_args()
model_fname = vars(args)['model_file']

# parse arguments for global variables 
SPECIES_ID = model_fname.split('/')[-1] # ID is model filename minus directory
SPECIES_ID = SPECIES_ID.split('.')[0] # get rid of extension
SPECIES_ID = SPECIES_ID.split('denovo_')[1]

day = datetime.now().strftime('%d_%m_%Y')
logging.basicConfig(filename='step4_{}_{}.log'.format(SPECIES_ID,day), level=logging.INFO, filemode='w')
logger = logging.getLogger(__name__)

logger.info('BEGIN STEP 4')

# modified for Rivanna: read in the models
model_dict = {}
model_dict[SPECIES_ID] = cobra.io.load_json_model(model_fname)

# load curated models
os.chdir(model_path)
chominis = cobra.io.read_sbml_model('chominis_vanee2010.xml')
tgondii = cobra.io.read_sbml_model('tg_tymoshenko2015.xml')
leish = cobra.io.read_sbml_model('iAC560_leishmania_Chavali2008.xml')
os.chdir(data_path)
pf_curated = cobra.io.load_json_model("iPfal18.json")

logger.info('loaded models')

# extract biomass reactions
pf_bm_mets = [x.copy() for x in pf_curated.reactions.get_by_id('biomass').metabolites]
ch_bm_mets = [x.copy() for x in chominis.reactions.get_by_id('R_biomass_target').metabolites]
leish_bm_mets = [x.copy() for x in leish.reactions.get_by_id('BIOMASS_LM3').metabolites]
tg_bm_mets = [x.copy() for x in tgondii.reactions.get_by_id('Biomass').metabolites]
tgondii_also_add = [tgondii.reactions.get_by_id('DNA').copy(),
					tgondii.reactions.get_by_id('RNA').copy(),
                    tgondii.reactions.get_by_id('GAM').copy(),
                    tgondii.reactions.get_by_id('Acyl').copy(),
                    tgondii.reactions.get_by_id('Proteins').copy(),
                    tgondii.reactions.get_by_id('Lipids').copy(),
                    tgondii.reactions.get_by_id('EMs').copy()]
tgondii_also_add_mets = [x.copy() for x in [y.reactants for y in tgondii_also_add]]
tgondii_also_add_mets = [item for sublist in tgondii_also_add_mets for item in sublist]

# get biomass precursor names
pf_biomass_mets_names = [x.name for x in pf_bm_mets]
chominis_biomass_mets_names = [x.name for x in ch_bm_mets]
leish_biomass_mets_names = [x.name for x in leish_bm_mets]
tgondii_also_add_mets_names = [x.name for x in tgondii_also_add_mets]

# make table of biomass precursor names
max_length = max(len(pf_biomass_mets_names),len(chominis_biomass_mets_names),
				 len(leish_biomass_mets_names),len(tgondii_also_add_mets_names))
pf_biomass_mets_ids = sorted(pf_biomass_mets_names, key=str.lower) + ['' for _ in \
range(max_length-len(pf_biomass_mets_names))]
for index, x in enumerate(pf_biomass_mets_ids):
    if '_' in x: pf_biomass_mets_ids[index] = x.replace('_','-')

chominis_biomass_mets_ids = sorted(chominis_biomass_mets_names, key=str.lower) + \
['' for _ in range(max_length-len(chominis_biomass_mets_names))]
leish_biomass_mets_ids = sorted(leish_biomass_mets_names, key=str.lower) + \
['' for _ in range(max_length-len(leish_biomass_mets_names))]
for index, x in enumerate(leish_biomass_mets_ids):
    if x.startswith('M_'):
        leish_biomass_mets_ids[index] = x[2:]
        x = leish_biomass_mets_ids[index]
    if x.endswith('_'):
        leish_biomass_mets_ids[index] = x[:-1]
        x = leish_biomass_mets_ids[index]
    if '__' in x:
        leish_biomass_mets_ids[index] = x.split('__')[0]
        x = leish_biomass_mets_ids[index]
    if '_' in x and x.split('_')[0] == x.split('_')[1]:
        leish_biomass_mets_ids[index] = x.split('_')[0]
        x = leish_biomass_mets_ids[index]
    if 'Cysteine' not in x and '_C' in x:
        leish_biomass_mets_ids[index] = x.split('_C')[0]
        x = leish_biomass_mets_ids[index]
    if '_' in x: leish_biomass_mets_ids[index] = x.replace('_','-')

tgondii_biomass_mets_ids = sorted(tgondii_also_add_mets_names, key=str.lower) + \
['' for _ in range(max_length-len(tgondii_also_add_mets_names))]
df = pd.DataFrame( {'iPfal17': pf_biomass_mets_ids,
					'C. hominis 2010': chominis_biomass_mets_ids,
                  	'T. gondii 2015': tgondii_biomass_mets_ids, 
                  	'Leishmania 2008': leish_biomass_mets_ids })
df.to_csv('biomass_components_{}.csv'.format(str(date.today())), header=True, index=True)

logger.info('wrote biomass components')

# get biomass precursor ids
pf_biomass_mets_ids = [x.id for x in pf_bm_mets]
chominis_biomass_mets_ids = [x.id for x in ch_bm_mets]
leish_biomass_mets_ids = [x.id for x in leish_bm_mets]
tgondii_also_add_mets_ids = [x.id for x in tgondii_also_add_mets]
lower_pf = set([met.lower().split('_')[0] for met in pf_biomass_mets_ids])
lower_lm = set([met.lower().split('_')[0] for met in leish_biomass_mets_ids])
lower_ch = set([met.lower().split('_')[0] for met in chominis_biomass_mets_ids])
biomass_mets = list(lower_ch.intersection(lower_lm).intersection(lower_pf))
biomass_mets.append('biomass')
biomass_mets.append('protein')
biomass_mets.append('lipid') # plus DAG begause in all reactions

logger.info('made biomass mets')

for species, model in model_dict.items():
    #biomass
    new_rxn = Reaction()
    met_dict = dict()
    for met in pf_curated.reactions.biomass.metabolites:
        if met.id.lower().split('_')[0] in biomass_mets:
            met_dict[met] = pf_curated.reactions.biomass.metabolites[met]
    new_rxn.add_metabolites(met_dict)
    new_rxn.name = 'generic_biomass'
    new_rxn.id = 'generic_biomass'
    new_rxn.lower_bound = 0.
    new_rxn.upper_bound = 1000.
    model.add_reactions([new_rxn])
    
    logger.info('added biomass')
    
    # protein production, transport, exchange
    # biomass 'sink', lipid exchange
    for x in ['Protein', 'Protein_t','DM_biomass_c','EX_lipid_c','Lipid_prod']: #'Protein_ex',
        new_rxn = pf_curated.reactions.get_by_id(x).copy()
        for met in new_rxn.metabolites:
            if met.id not in [x.id for x in model.metabolites]:
                model.add_metabolites([met])
        model.add_reactions([new_rxn])
        if len(new_rxn.genes) > 0:
            for gene in new_rxn.genes:
                model.genes.remove(gene)
    
    # lipid production
    for r in pf_curated.metabolites.get_by_id('all_dgl_c').reactions:
        if r.id != 'Lipid_prod':
            new_rxn = pf_curated.reactions.get_by_id(r.id).copy()
            for met in new_rxn.metabolites:
                if met.id not in [x.id for x in model.metabolites]:
                    model.add_metabolites([met])
            model.add_reactions([new_rxn])
            if len(new_rxn.genes) > 0:
                for gene in new_rxn.genes:
                    model.genes.remove(gene)
    new_rxn = Reaction()
    met_dict = {pf_curated.metabolites.get_by_id('all_dgl_c'): -1,
        pf_curated.metabolites.get_by_id('lipid_c'): +1}
    new_rxn.add_metabolites(met_dict)
    new_rxn.name = 'Lipid_prod'
    new_rxn.id = 'Lipid_prod'
    new_rxn.lower_bound = 0.
    new_rxn.upper_bound = 1000.
    for met in new_rxn.metabolites:
        if met.id not in [x.id for x in model.metabolites]:
            model.add_metabolites([met])
    model.add_reactions([new_rxn])
    
    logger.info('added accessory reactions')

    model.objective = model.reactions.get_by_id('generic_biomass')
    
    logger.info('set objective')

    # add exchange reactions for all extracellular metabolites to accelerate gapfilling
    for met in model.metabolites:
        if met.id.endswith('_e'):
            if 'EX_'+met.id not in model.reactions:
                model.add_boundary(met, type = 'exchange')
    for rxn in model.reactions:
        if rxn.id.startswith('EX_'):
            rxn.lower_bound = -1000.
            rxn.upper_bound = 1000.

    logger.info('fixed exchange')
    
    model_dict[species] = model
    ### SAVE MODEL
    os.chdir(model_path)

    if 'hb_c' in [m.id for m in model.metabolites]:
        logger.info('HEMOGLOBIN PRESENT')
    os.chdir(data_path)
    cobra.io.save_json_model(model, "./with_biomass_"+species+".json")


