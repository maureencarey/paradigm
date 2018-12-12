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

logging.info('BEGIN STEP 3 for species:')
logging.info(SPECIES_ID)

species = SPECIES_ID
x = species
logging.info('x = ')
if '.tsv' in x:
    logging.info('.tsv in the annotations file string, this might cause problems')
    species = x.split('.tsv')[0]
else:
    species = x
if '_BiGG' in x:
    logging.info('_BiGG in the annotations file string, this might cause problems')
    species = species.split('_BiGG')[0]
logging.info(x)

logging.info("DIY2_"+species+".json")
logging.info('______ FINAL OUTPUT OF THIS FILE IS_____')
logging.info("final_denovo_"+species+".json")
