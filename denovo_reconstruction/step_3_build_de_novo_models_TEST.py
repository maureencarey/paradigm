import cobra
import pandas as pd
import os
import subprocess
import glob
import json
from cobra import Model, Reaction, Metabolite
import helper_functions_1 as hf
import argparse
import logging
from datetime import datetime

parser = argparse.ArgumentParser(description='Read in the species annotation filename')
parser.add_argument('annotation_file')
args = parser.parse_args()
annotation_fname = vars(args)['annotation_file']

# parse arguments for global variables
SPECIES_ID = annotation_fname.split('/')[-1] # ID is annotation filename minus directory
SPECIES_ID = SPECIES_ID.split('_BiGG.')[0] # get rid of extension

logging.basicConfig(filename='log_file.log', level=logging.INFO)
logger = logging.getLogger(__name__)

logging.info('BEGIN STEP 3')
species = SPECIES_ID
logger.info('SPECIES ID is'+SPECIES_ID)
logger.info('____Input MODEL NAME IS ____'+annotation_fname)
logger.info('____OUTPUT MODEL NAME IS final_denovo_'+species+'.json____')
