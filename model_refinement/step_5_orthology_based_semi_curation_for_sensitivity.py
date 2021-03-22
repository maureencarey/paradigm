## Input = with_biomass_denovo_SPECIES.json
## ONLY DO PLSASMODIUMS

import cobra
import pandas as pd
import os
import subprocess
import glob
import json
from cobra import Model, Reaction, Metabolite
#import helper_functions_2 as hf
import sys
sys.path.append(os.path.abspath("/home/mac9jc/paradigm/"))
import helper_functions as hf
import argparse
import logging
from datetime import datetime
from datetime import date


data_path = "/home/mac9jc/paradigm/data"
model_path = "/home/mac9jc/paradigm/models"
log_file_path = "/home/mac9jc/paradigm/model_generation_logs"
os.chdir(log_file_path)

parser = argparse.ArgumentParser(description='Read in the species model')
parser.add_argument('model_file')
args = parser.parse_args()
model_fname = vars(args)['model_file']

# parse arguments for global variables 
SPECIES_ID = model_fname.split('/')[-1] # ID is model filename minus directory
SPECIES_ID = SPECIES_ID.split('.')[0] # get rid of extension

day = datetime.now().strftime('%d_%m_%Y')
logging.basicConfig(filename='step5_for_sensitivity_{}_{}.log'.format(SPECIES_ID,day), level=logging.INFO, filemode='w')
logger = logging.getLogger(__name__)

SPECIES_ID = SPECIES_ID.split('with_biomass_')[1]
logger.info('BEGIN STEP 5')

os.chdir(model_path)
iPfal19 = cobra.io.load_json_model("iPfal21.json")

logger.info(SPECIES_ID)
# modified for Rivanna: read in the models
model = cobra.io.load_json_model(model_fname)
logger.info('loaded model')

os.chdir(data_path)
mapping = pd.read_csv("plasmodium_orthology_conversion_release44.tsv", sep='\t')
mapping['Organism'] = mapping['Organism'].str.replace(".", "")
mapping['Organism'] = mapping['Organism'].str.replace("strain", "")
mapping['Organism'] = mapping['Organism'].str.replace("Strain", "")
mapping['Organism'] = mapping['Organism'].str.replace(" ", "")
mapping['Organism'] = mapping['Organism'].str.replace("Sal-1", "Sal1")
mapping['Organism'] = mapping['Organism'].str.replace("nilgiri","Nilgiri")
if SPECIES_ID not in mapping['Organism'].unique():
    logger.info(SPECIES_ID+'not in orthology mapping dataframe, aka we are not converting via orthology (printing infeasible to catch this error)')
else:
    logger.info('loaded mapping')
add_these = ['pe_prod1', 'pe_prod10', 'pe_prod11', 'pe_prod12', 'pe_prod13', 'pe_prod14',
             'pe_prod15', 'pe_prod16', 'pe_prod17', 'pe_prod18', 'pe_prod19', 'pe_prod20',
             'pe_prod21', 'pe_prod22','pe_prod2','pe_prod3','pe_prod4','pe_prod5',
             'pe_prod6','pe_prod7','pe_prod8', 'pe_prod9', 'pg_prod1', 'pg_prod2', 
             'pg_prod3','pg_prod4', 'pg_prod5', 'pg_prod6', 'pg_prod7','pi_prod1',
             'pi_prod2','pi_prod3','pi_prod4',#'progly_t',
             'ps_prod1', 'ps_prod2', 'ps_prod3','ps_prod4', 'ps_prod5','ps_prod6',
             'ps_prod7','lipid1', 'lipid10', 'lipid2','lipid3','lipid4', 'lipid5',
             'lipid6', 'lipid7', 'lipid8','lipid9','PA120tap','PA140tap', 'PA141tap',
             'PA160tap', 'PA161tap', 'PA180tap', 'PA181tap','2agpe120t',
             '2agpe140t','2agpe141t','2agpe160t','2agpe161t', '2agpe180t',
             '2agpe181t', '2agpg120t','2agpg140t','2agpg141t','2agpg160t', 
             '2agpg161t','2agpg180t', '2ddecg3pt','acoa_prod1','acoa_prod10', 
             'acoa_prod11', 'acoa_prod12', 'acoa_prod13', 'acoa_prod14',
             'acoa_prod15','acoa_prod16','acoa_prod2', 'acoa_prod3', 'acoa_prod4',
             'acoa_prod5', 'acoa_prod6', 'acoa_prod7', 'acoa_prod8', 'acoa_prod9',
             'acylpg_prod1', 'acylpg_prod2','acylpg_prod3','acylpg_prod4','acylpg_prod5',
             'acylpg_prod6','acylpg_prod7',  'dolichol_t',#'TTDCAt',#'FEROpp',
             'XOLEST2te','PCt','pail_uptake','pc_prod', 'Htap','Htmt','2agpg181t',
             'CHSTEROLt','HBtr','SM_host','HMBZex','Hfv','O2St','O2Stm']
# bm models already have biomass biomass_generic 'protein_t_e', 'dgl_prod1','dgl_prod2',
# 'Lipid_prod','dgl_prod3', 'dgl_prod4', 'dgl_prod5','dgl_prod6', 'dgl_prod7','Protein',

columns = ['species','starting_genes','reactions_added', 'mets_added','genes_added']
modifications_ortho = pd.DataFrame(index = [0], columns=columns)

logger.info('beginning model loop')

# keep bounds from universal model

x1 = len(model.reactions)
y1 = len(model.genes)
z1 = len(model.metabolites)

# add reactions that are necessary for biomass (aggregation and transport rxns)
for x in add_these:
    rxn = iPfal19.reactions.get_by_id(x).copy()
    if rxn.gene_reaction_rule:
        logger.info(rxn.id)
        logger.info('has the following GPR:')
        logger.info(rxn.gene_reaction_rule)
    for met in rxn.metabolites:
        if met.id not in [m.id for m in model.metabolites]:
            model.add_metabolites([met.copy()])
    model.add_reactions([rxn])

logger.info('added basic reactions')

x2 = len(model.reactions)
y2 = len(model.genes)
z2 = len(model.metabolites)

species_specific_mapping = mapping[mapping['Organism'] == SPECIES_ID].copy()
logger.info('pruned mapping file')

for index, row in species_specific_mapping.iterrows():
    new_gene = row['Gene ID']
    if ',' not in row['Input Ortholog(s)']:
        gene = row['Input Ortholog(s)'].strip()
        if gene in [x.id for x in iPfal19.genes]:
            gene = iPfal19.genes.get_by_id(gene)
            rxn_to_add = gene.reactions
            for rxn in rxn_to_add:
                if rxn.id not in [x.id for x in model.reactions]:
                    rxn2 = rxn.copy()
                    mets = [x.metabolites for x in [rxn2]]
                    all_keys = set().union(*(d.keys() for d in mets))
                    for key in all_keys:
                        if key.id not in [m.id for m in model.metabolites]:
                            model.add_metabolites([key.copy()])
                    model.add_reactions([rxn2])
                    model.reactions.get_by_id(rxn2.id).gene_reaction_rule = new_gene
                else: # reaction is in model
                    if new_gene not in model.reactions.get_by_id(rxn.id).gene_reaction_rule:
                        if model.reactions.get_by_id(rxn.id).gene_reaction_rule == '':
                            model.reactions.get_by_id(rxn.id).gene_reaction_rule = new_gene
                        else: model.reactions.get_by_id(rxn.id).gene_reaction_rule = \
                            model.reactions.get_by_id(rxn.id).gene_reaction_rule + ' or ' +new_gene
#                     else:logger.info('gene already there')
    else:
        for i in range(0,len(row['Input Ortholog(s)'].split(', ')),1):
            gene = row['Input Ortholog(s)'].split(', ')[i].strip()
            if gene in [x.id for x in iPfal19.genes]:
                gene = iPfal19.genes.get_by_id(gene)
                rxn_to_add = gene.reactions
                for rxn in rxn_to_add:
                    if rxn.id not in [x.id for x in model.reactions]:
                        rxn2 = rxn.copy()
                        mets = [x.metabolites for x in [rxn2]]
                        all_keys = set().union(*(d.keys() for d in mets))
                        for key in all_keys:
                            if key.id not in [m.id for m in model.metabolites]:
       	                        model.add_metabolites([key.copy()])
                        model.add_reactions([rxn2])
                        model.reactions.get_by_id(rxn2.id).gene_reaction_rule = new_gene
                    else:
                        if new_gene not in model.reactions.get_by_id(rxn.id).gene_reaction_rule:
                            if model.reactions.get_by_id(rxn.id).gene_reaction_rule == '':
                                model.reactions.get_by_id(rxn.id).gene_reaction_rule = new_gene
                            else: model.reactions.get_by_id(rxn.id).gene_reaction_rule = \
                                model.reactions.get_by_id(rxn.id).gene_reaction_rule + ' or ' +new_gene
#                             else: #gene already there
if len(model.reactions) != len(set(model.reactions)):
    logger.info('duplicate reactions present')

model.repair()

if SPECIES_ID.startswith('Pfalciparum3D7'):
    s = 1
else:
    t = len(model.genes)
    gene_list = [x.id for x in model.genes if x.id.startswith('PF3D7')]
    gene_list_genes = [x for x in model.genes if x.id.startswith('PF3D7')]
    if len(gene_list) > 0:
        logger.info(len(gene_list))
        for x in gene_list_genes:
            if len(x.reactions) == 0:
                model.genes.remove(x)
            else:
                logging.info(x.id+' remains in model and associated with '+[r.id for r in x.reactions])
                cobra.manipulation.delete.remove_genes(model, gene_list, remove_reactions=False)
        if (t==len(model.genes)) or len([x.id for x in model.genes if x.id.startswith('PF3D7')]) > 0:
            logger.info('ERROR: (STEP 5) DELETE PF3D7 GENES STEP DID NOT WORK')
    gene_list = [x.id for x in model.genes if x.id in ['mal_mito_1','mal_mito_2','mal_mito_3']]
    if len(gene_list) > 0:
        cobra.manipulation.delete.remove_genes(model, gene_list, remove_reactions=False)
        if len([x.id for x in model.genes if x.id in ['mal_mito_1','mal_mito_2','mal_mito_3']]) > 0:
            logger.info('ERROR: (STEP 5) DELETE PF3D7 mal_mito GENES STEP DID NOT WORK')

# save the number of modifications
x3 = len(model.reactions)
y3 = len(model.genes)
z3 = len(model.metabolites)
row_index = 0#modifications_ortho.species == SPECIES_ID
modifications_ortho.loc[row_index,'species'] = SPECIES_ID
modifications_ortho.loc[row_index,'starting_genes'] = y1
modifications_ortho.loc[row_index,'reactions_added'] = len(model.reactions) - x1 
modifications_ortho.loc[row_index,'mets_added'] = len(model.metabolites) - z1
modifications_ortho.loc[row_index,'genes_added'] = len(model.genes) - y1 

logger.info('headed into exchanges')
# add exchanges in preparation for gapfilling
l = list()
for rxn in model.reactions:
    for met in rxn.metabolites:
        if met.id.endswith('_e'):
            l.append(met)
for met in l:
    if 'EX_'+met.id not in [r.id for r in model.reactions]:
        model.add_boundary(met, type = "exchange")
model.objective = 'biomass' # NOT USING GENERIC BIOMASS NOW

model.repair()
os.chdir(model_path)
cobra.io.save_json_model(model, "for_sensitivity_ortho_"+SPECIES_ID+".json")
cobra.io.write_sbml_model(model, "for_sensitivity_ortho_"+SPECIES_ID+".xml")

if 'hb_c' in [m.id for m in model.metabolites]:
    logger.info('HEMOGLOBIN PRESENT')

logger.info('saved model')
logger.info('no genes should print here (if a number, then your gene_IDs were not updated):')
if SPECIES_ID != 'Pfalciparum3D7':
    for x in model.genes:
        if x.id.startswith('PF3D7'):
            logger.info(x.id)
    logger.info('k, moving on')

# SKIP FOR NOW: FIND/ REMOVE DUPLICATE REACTIONS
logger.info('duplicate reactions by formula')
duplicates = dict()
temp_dict = dict() # get products and reactants for every reaction
for rxn in model.reactions:
    rxn_dict = dict()
    check_rxn_products = rxn.products
    check_rxn_reactants = rxn.reactants
    rxn_dict['reactants'] = [x.id for x in check_rxn_reactants]
    rxn_dict['products'] = [x.id for x in check_rxn_products]
    temp_dict[rxn.id] = rxn_dict
    
    for key in temp_dict.keys():
        if key != rxn.id:
            if temp_dict[rxn.id]['reactants'] == temp_dict[key]['reactants'] and \
            temp_dict[rxn.id]['products'] == temp_dict[key]['products']:
                if rxn.id not in duplicates.keys():
                    duplicates[rxn.id] = key
                elif duplicates[rxn.id] == key or key in duplicates[rxn.id]:
                    continue
                else:
                    duplicates[rxn.id] = duplicates[rxn.id]+', '+key
logger.info(duplicates)

model.repair()

os.chdir(data_path)
modifications_ortho.to_csv('./orthology_modifications_plasmodium_for_sensitivity_{0}.csv'.format(SPECIES_ID))
