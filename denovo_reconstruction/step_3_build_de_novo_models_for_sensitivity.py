import cobra
import pandas as pd
import os
import subprocess
import glob
import json
from cobra import Model, Reaction, Metabolite
#import helper_functions_1 as hf
import sys
sys.path.append(os.path.abspath("/home/mac9jc/paradigm/"))
import helper_functions as hf
import argparse
import logging
from datetime import datetime


########
# this script is different from step_3_build_de_novo_models.py because  'bad' reactions, 
# or reactions in an unaccepted compartment, are removed and NOT switched into the cytosol
# also only planning on running this for T. gondii, P. falciparum, and P. berghei since they
# have genome-wide essentiality screens
########


parser = argparse.ArgumentParser(description='Read in the species annotation filename')
parser.add_argument('annotation_file')
args = parser.parse_args()
annotation_fname = vars(args)['annotation_file']

# parse arguments for global variables
SPECIES_ID = annotation_fname.split('/')[-1] # ID is annotation filename minus directory
SPECIES_ID = SPECIES_ID.split('_BiGG.')[0] # get rid of extension

# set up log file
log_file_path =	"/home/mac9jc/paradigm/model_generation_logs"
os.chdir(log_file_path)
day = datetime.now().strftime('%d_%m_%Y')
logging.basicConfig(filename='step3_for_sensitivity_{}_{}.log'.format(SPECIES_ID,day), level=logging.INFO, filemode='w')
logger = logging.getLogger(__name__)

# begin
logging.info('BEGIN STEP 3')
logging.info(SPECIES_ID)
data_path = "/home/mac9jc/paradigm/data"
model_path = "/home/mac9jc/paradigm/models"
og_path = "/home/mac9jc/paradigm/"

# universal reaction bag for model generation
os.chdir(model_path)
universal_model = cobra.io.read_sbml_model('universal/universal_model_updated.xml')
universal_model = hf.update_universal_model(universal_model)
os.chdir(data_path)
logging.info('loaded universal')

# compartments in universal reaction bag
compartment_options = list()
for string in [met.id for met in universal_model.metabolites]:
    if string.startswith('EX_') or string.startswith('DM_') or string.startswith('SK_'):
        s = 'asedf'
    else:
        compartment_options.append('_'+string.split('_')[len(string.split('_'))-1])
compartment_options = set(compartment_options)
compartment_options = [x.split()[0] for x in compartment_options]
logging.info(compartment_options)

# get metabolites involved in each reaction in universal model in dictionary format
universal_dict = dict() 
for rxn in universal_model.reactions:
    rxn_dict = dict()

    check_rxn_products = rxn.products
    check_rxn_reactants = rxn.reactants
    all_mets = rxn.metabolites
    compart = set([hf.get_comp(universal_model,x.id) for x in all_mets])

    rxn_dict['reactants'] = [hf.met_ids_without_comp(universal_model,x.id) for x in check_rxn_reactants]
    rxn_dict['products'] = [hf.met_ids_without_comp(universal_model,x.id) for x in check_rxn_products]
    rxn_dict['compartment'] = list(compart)

    universal_dict[rxn.id] = rxn_dict
    # universal_dict = all universal model reactions, mapped to a dictionary
    # containing its compartment, products and reactants

# map reactions to duplicate reactions in different compartments
universal_dict_with_alts = universal_dict.copy() 
for reaction, data in universal_dict_with_alts.items():
    alternative_reactions = dict()
    data_with_options = data.copy()
    for potential_rxn,potential_data in universal_dict.items():
        if potential_rxn != reaction:
            if potential_data['reactants'] == data['reactants'] and \
            potential_data['products'] == data['products']:
                alternative_reactions[potential_rxn] = potential_data['compartment']

    data_with_options['alternative_reactions'] = alternative_reactions
    universal_dict_with_alts[reaction] = data_with_options

# database mapping - this is for release 42, must update if using a different EuPathDB release
plasmodb = ["PadleriG01","PbergheiANKA","PbillcollinsiG01","PblacklockiG01","Pchabaudichabaudi","PcoatneyiHackeri","PcynomolgiB","PcynomolgiM","Pfalciparum3D7","Pfalciparum7G8","PfalciparumCD01","PfalciparumDd2","PfalciparumGA01","PfalciparumGB4","PfalciparumGN01","PfalciparumHB3","PfalciparumIT","PfalciparumKE01","PfalciparumKH01","PfalciparumKH02","PfalciparumML01","PfalciparumSD01","PfalciparumSN01","PfalciparumTG01","PfragileNilgiri","PgaboniG01","PgaboniSY75","Pgallinaceum8A","PinuiSanAntonio1","PknowlesiH","PknowlesiMalayanPk1A","PmalariaeUG01","PovalecurtisiGH01","PpraefalciparumG01","PreichenowiCDC","PreichenowiG01","PrelictumSGS1-like","PvinckeipetteriCR","Pvinckeivinckeivinckei","Pvivax-likePvl01","PvivaxP01","PvivaxSal1","Pyoeliiyoelii17X","Pyoeliiyoelii17XNL","PyoeliiyoeliiYM"]
toxodb = ["CcayetanensisCHN_HEN01", "CsuisWienI","EacervulinaHoughton", "EbrunettiHoughton", "EfalciformisBayerHaberkorn1970", "EmaximaWeybridge", "EmitisHoughton", "EnecatrixHoughton", "EpraecoxHoughton", "EtenellaHoughton", "HhammondiHH34", "NcaninumLIV", "SneuronaSN3", "SneuronaSOSN1", "TgondiiARI", "TgondiiFOU", "TgondiiGAB2-2007-GAL-DOM2", "TgondiiGT1", "TgondiiMAS", "TgondiiME49", "Tgondiip89", "TgondiiRH", "TgondiiRUB", "TgondiiTgCatPRC2", "TgondiiVAND", "TgondiiVEG"]

# load annotation data file
os.chdir(data_path)
os.chdir("./diamond_output_BiGG")
annotations_file = pd.read_table(annotation_fname)
annotations_file.columns = ['query_gene', 'BiGG_gene', 'pident', 'length', 'mismatch', 'gapopen','qstart', 'qend', 'sstart', 'send', 'evalue', 'score']

# map annotations to BiGG functions for model generation
os.chdir(data_path)
gprs = pd.read_csv('./bigg_gprs.csv') # From CarveMe
gprs.reaction = [x[2:] for x in gprs.reaction]
gprs = gprs[gprs.reaction.isin([rxn.id for rxn in universal_model.reactions])] # updated from CarveMe
logging.info("loaded gprs")

# score reactions per annotations
scores_file = hf.reaction_scoring(annotations_file, gprs)
# carveme will maximize positive scores and minimize negative scores while maintaining a functional network
# we will just add any positive scoring reaction above a certain threshold and make it functional later
keep_scores = scores_file.loc[scores_file.score>10]
logging.info("done with scoring")

# make model from universal model
new_model = universal_model.copy()
new_model.name = SPECIES_ID
new_model.id = SPECIES_ID
starting_rxn_ct = len(new_model.reactions)

if len(keep_scores.reaction) == len(set(keep_scores.reaction)):
    
    rxns_to_add = dict()
    scores_for_rxns = dict()
    for index, row in keep_scores.iterrows():
        rxns_to_add[row['reaction']] = row['GPR']
        scores_for_rxns[row['reaction']] = row['score']
    new_model.remove_reactions([rxn for rxn in new_model.reactions if rxn.id not in rxns_to_add.keys()])
    new_model.repair()

    if not [rxn.id for rxn in new_model.reactions if rxn.gene_reaction_rule != '']:
        for rxn in new_model.reactions:
                if rxn.gene_reaction_rule == '':
                    new_model.reactions.get_by_id(rxn.id).gene_reaction_rule = rxns_to_add[rxn.id]
                    new_model.reactions.get_by_id(rxn.id).notes['CarveMe score'] = {rxns_to_add[rxn.id] : scores_file.loc[scores_file.reaction == rxn.id]['score'].values[0]}
    else:
        logging.info('error - some reactions already have GPRs')
    logging.info('made new model from universal')
    new_model.repair()

    if len(rxns_to_add.keys()) == len(new_model.reactions):
        if starting_rxn_ct <= len(rxns_to_add.keys()):
            logging.info('error with original model, reactions already removed')
    else:
        logging.info('error with reaction removal, resultant len(model.reactions) != rxns_to_keep')
else:
    logging.info('duplicate keep_scores.reaction')

if len(rxns_to_add.keys()) != len(new_model.reactions):
    logging.info('error in universal reaction pruning')

logging.info('made first draft model, with this many reactions:')
logging.info(len(new_model.reactions))
    
# remove duplciate reactions in mulitple compartments
total_compartments = ["_c","_e","_m","_ap","_fv","_k","_glc","_pm"]
# cytosol, extracellular, mitochondrdia, apicoplast, food vacuole, kinetoplast, glycosome, pseudomitochondria

if SPECIES_ID in plasmodb: # Plasmodium = cytosol, extracellular, mitochondrdia, apicoplast, food vacuole
    compartment = ["_c","_e","_m","_ap","_fv"]
elif SPECIES_ID in toxodb: # Toxoplasma = cytosol, extracellular, mitochondrdia, apicoplast
    compartment = ["_c","_e","_ap","_m"]
else:
    compartment = ["_c","_e"]

model = new_model
logging.info('finding good or bad reactions')
not_compartments = list(set(compartment_options) - set(compartment))

# get reactions that use/make at least one metabolite that is in an inappropariate compartment
good_rxns = list()
bad_rxns = list()
# not_compartments = [x+' ' for x in not_compartments]
for rxn_object in model.reactions: # if a reaction does not contain any bad compartments
    rxn_bad_counter = 0
    for x in not_compartments:
        if x in rxn_object.reaction or rxn_object.reaction.endswith(x):
            rxn_bad_counter = rxn_bad_counter + 1
    if rxn_bad_counter == 0:
        good_rxns.append(rxn_object.id)
    else:
        bad_rxns.append(rxn_object.id)

# how many bad rxns have genes associated with them that will be removed in this process?
rxn_gene_dict = dict()
for r_id in bad_rxns:
    rxn = model.reactions.get_by_id(r_id)
    i = 0
    if len(rxn.genes) > 0:
        for gene in rxn.genes:
            if len(gene.reactions) > 1: i = i + 1
    if i > 1: answer = 'yes' 
    else: answer = 'no'
    rxn_gene_dict[r_id] = {'genes':rxn.genes, 'gene_assoc_w_other_rxns?':answer}        
df_temp = pd.DataFrame.from_dict(rxn_gene_dict, orient='index')
df_temp.to_csv('/home/mac9jc/paradigm/data/results/genes_that_require_compartmentalization_'+SPECIES_ID+day+'.csv') 

# remove reactions
x = len(model.reactions)
y = len(model.metabolites)
model.remove_reactions(bad_rxns)
for gene in model.genes:
    if len(gene.reactions) == 0:
        cobra.manipulation.remove_genes(model,[gene.id])
model.repair()

if len(bad_rxns) != (x - len(model.reactions)):
    logging.info('error - reaction not removed properly')
x1 = x - len(model.reactions)
y1 = y - len(model.metabolites)
x1_2 = len(model.reactions)

# make sure all reactions can carry flux, problem with some versions of the universal model
for rxn in model.reactions:
    if rxn.lower_bound == 0 and rxn.upper_bound == 0:
        logging.info(rxn.id + ' has bounds == 0 in '+key)
        rxn.lower_bound = -1000.
        rxn.upper_bound = 1000.
        # NOTHING SHOULD PRINT - this was a problem in CarveMe
model, unused = hf.prune_unused_metabolites2(model)

# ADD IN THIS WARNING:
for c in not_compartments:
    if c in ([hf.get_comp(model,x.id) for x in model.metabolites]):
        logging.info('ERROR, UNACCEPTABLE COMPARTMENTS:')
        logging.info(c)

l2 = list()
for rxn in model.reactions:
    for suffix in [hf.get_comp(model,m.id) for m in rxn.metabolites]:
        l2.append(suffix)
logging.info('compartments:')
logging.info(set(l2))

model, unused = hf.prune_unused_metabolites2(model)
model.solver = 'glpk'

# check compartments
list_om= list()
list_om2= list()
for rxn in model.reactions:
    for m in rxn.metabolites:
        list_om.append(hf.get_comp(model,m.id))
for m in model.metabolites:
    list_om2.append(hf.get_comp(model,m.id))
if set(list_om) != set(list_om2):
    logging.info('error - extra compartments are present, pruning of unused metabolites did not work')

compartment_dict = {'c': 'cytoplasm', 'e': 'extracellular', 'm': 'mitochondrion', 'fv': 'food vacuole', 'ap':'apicoplast','k':'kinetoplast','glc':'glycosome', 'p':'periplasm','x':'peroxisome/glyoxysome','r':'endoplasmic reticulum','v':'vacuole','n':'nucleus','g':'golgi apparatus','u':'thylakoid','l':'lysosome','h':'chloroplast','f':'flagellum','s':'eyespot','im':'intermembrane space of mitochondria','cx':'carboxyzome','um':'thylakoid membrane','cm':'cytosolic membrane','i':'not_provided_by_bigg'}
for met in model.metabolites:
    if met.compartment == '':
        comp = hf.get_comp(model,met.id)
        comp = comp[1:]
        if comp in compartment_dict:
            comp_use = compartment_dict[comp]
        else: comp_use = 'no compartment'
        model.metabolites.get_by_id(met.id).compartment = comp_use

for gene in model.genes:
    if len(gene.reactions) == 0:
        cobra.manipulation.remove_genes(model,[gene.id])
model.repair()
os.chdir(model_path)
cobra.io.save_json_model(model, "for_sensitivity_denovo_"+SPECIES_ID+".json")
