import cobra
import os
import pandas as pd
from cobra.core import Gene, Metabolite, Reaction
import requests
import time
import logging
#import helper_functions_3 as hf3
import sys
sys.path.append(os.path.abspath("/home/mac9jc/paradigm/"))
import helper_functions as hf

model_path = "/home/mac9jc/paradigm/models/"
log_path = "/home/mac9jc/paradigm/model_generation_logs/"

os.chdir(log_path)
logging.basicConfig(filename='update_universal_round1b.log', level=logging.INFO, filemode='w')
logger = logging.getLogger(__name__)

os.chdir(model_path)
model = cobra.io.load_json_model('universal_model_updated.json')

compartment_dict = {'c': 'cytoplasm', 'e': 'extracellular', 'm': 'mitochondrion', 'fv': 'food vacuole', 'ap':'apicoplast','k':'kinetoplast','glc':'glycosome', 'p':'periplasm','x':'peroxisome/glyoxysome','r':'endoplasmic reticulum','v':'vacuole','n':'nucleus','g':'golgi apparatus','u':'thylakoid','l':'lysosome','h':'chloroplast','f':'flagellum','s':'eyespot','im':'intermembrane space of mitochondria','cx':'carboxyzome','um':'thylakoid membrane','cm':'cytosolic membrane','i':'not_provided_by_bigg'}

for met in model.metabolites:
    comp = hf.get_comp(model,met.id)    
    comp = comp[1:]
    if comp in compartment_dict:
        comp_use = compartment_dict[comp]
    else: comp_use = 'no compartment'
    model.metabolites.get_by_id(met.id).compartment = comp_use

os.chdir(model_path)
cobra.io.save_json_model(model,  'universal_model_updated2.json')
cobra.io.write_sbml_model(model,  'universal_model_updated2.xml')

