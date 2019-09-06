import cobra
import os
from cobra import Model, Reaction, Metabolite
import pandas as pd
import requests
import logging
from datetime import datetime
import sys
sys.path.append(os.path.abspath("/home/mac9jc/paradigm/"))
import helper_functions as hf

og_path = "/home/mac9jc/paradigm"
data_path = "/home/mac9jc/paradigm/data"
model_path = "/home/mac9jc/paradigm/models"
day = datetime.now().strftime('%d_%m_%Y')

os.chdir(model_path)
pf_model = cobra.io.read_sbml_model("iPfal19_without_annotation.xml")

pf_model.annotation["taxonomy"] = "36329"
pf_model.annotation["genome"] = "https://plasmodb.org/common/downloads/Current_Release/Pfalciparum3D7/fasta/data/PlasmoDB-44_Pfalciparum3D7_Genome.fasta"
pf_model.annotation["DOI"] = "pending"
pf_model.annotation["authors"] = 'Maureen Carey, mac9jc@virginia.edu'
# cobrapy cannot currently handle the following:
#pf_model.annotation["authors"] =  [
#{"familyName": "Carey","givenName": "Maureen","organisation": "University of Virginia","email": "mac9jc@virginia.edu"},
#{"familyName": "Untaroiu","givenName": "Ana","organisation": "University of Virginia","email": "amu4pv@virginia.edu"}, 
#{"familyName": "Plata","givenName": "German","organisation": "Columbia University","email": "gap2118@columbia.edu"}]
pf_model.annotation["species"] = "Plasmodium falciparum"
pf_model.annotation["strain"] = "3D7"
pf_model.annotation["tissue"] = "parasite in the asexual blood-stage"
pf_model.annotation["terms_of_distribution"] = "CC-BY" 
pf_model.annotation["created"] = day
pf_model.annotation["sbo"] = "SBO:0000624"
pf_model.annotation["curation"] = ['DOI: 10.1038/msb.2010.60', 'DOI: 10.1186/s12864-017-3905-1', 'DOI: 10.1186/s12859-019-2756-y', 'unpublished by Maureen Carey']
pf_model.annotation["genedb"] = "Pfalciparum"

cobra.io.write_sbml_model(pf_model, "iPfal19_with_annotation.xml")
