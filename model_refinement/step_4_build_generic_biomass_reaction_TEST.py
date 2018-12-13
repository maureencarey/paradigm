## Input = final_denovo_SPECIES.json

import argparse
import logging

parser = argparse.ArgumentParser(description='Read in the species model')
parser.add_argument('model_file')
args = parser.parse_args()
model_fname = vars(args)['model_file']

# parse arguments for global variables 
SPECIES_ID = model_fname.split('/')[-1] # ID is model filename minus directory
SPECIES_ID = SPECIES_ID.split('.')[0] # get rid of extension
SPECIES_ID = SPECIES_ID.split('denovo_')[1]

logging.basicConfig(filename='log_file.log', level=logging.INFO)
logger = logging.getLogger(__name__)

logger.info('BEGIN STEP 4')
species = SPECIES_ID
logger.info('SPECIES ID is'+SPECIES_ID)
logger.info('____Input MODEL NAME IS ____'+model_fname)
logger.info('____OUTPUT MODEL NAME IS with_biomass_'+species+'.json____')

