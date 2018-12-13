## Input = with_biomass_denovo_SPECIES.json
## ONLY DO PLSASMODIUMS

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
from datetime import datetime

parser = argparse.ArgumentParser(description='Read in the species model')
parser.add_argument('model_file')
args = parser.parse_args()
model_fname = vars(args)['model_file']

# parse arguments for global variables 
SPECIES_ID = model_fname.split('/')[-1] # ID is model filename minus directory
SPECIES_ID = SPECIES_ID.split('.')[0] # get rid of extension
SPECIES_ID = SPECIES_ID.split('with_biomass_')[1]

logging.basicConfig(filename='log_file.log', level=logging.INFO)
logger = logging.getLogger(__name__)
logger.info('BEGIN STEP 5')

logger.info(SPECIES_ID)
species = SPECIES_ID
logger.info('output model name is = ')
logger.info("ortho_"+species+".json")

