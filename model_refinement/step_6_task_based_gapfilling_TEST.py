# Input ortho_speciesname.json and iPfal18.json

import cobra
import pandas as pd
import os
import subprocess
import glob
import json
from cobra import Model, Reaction, Metabolite
import helper_functions_2 as hf
import argparse
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
SPECIES_ID_old = SPECIES_ID
if 'final_denovo_' in SPECIES_ID:
    SPECIES_ID = SPECIES_ID.split('final_denovo_')[1]
if 'with_biomass_' in SPECIES_ID:
    SPECIES_ID = SPECIES_ID.split('with_biomass_')[1]
if 'ortho_' in SPECIES_ID:
    SPECIES_ID = SPECIES_ID.split('ortho_')[1]
if 'DIY' in SPECIES_ID:
    SPECIES_ID = SPECIES_ID[5:]
if SPECIES_ID == 'iPfal18':
    SPECIES_ID = 'Pfalciparum3D7'

logging.basicConfig(filename='log_file.log', level=logging.INFO)
logger = logging.getLogger(__name__)
logger.info('BEGIN STEP 6')
logger.info('SPECIES ID is')
logger.info(SPECIES_ID)
    
species = SPECIES_ID
logger.info('FINAL MODEL NAME IS gf_'+SPECIES_ID_old+'.json')
    
    
