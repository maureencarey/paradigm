
import cobra
import os
import subprocess
import glob
import json
from cobra import Model, Reaction, Metabolite
import sys
sys.path.append(os.path.abspath("/home/mac9jc/paradigm/"))
import helper_functions as hf
import logging
from datetime import datetime

log_file_path =	"/home/mac9jc/paradigm/model_generation_logs"
model_path = "/home/mac9jc/paradigm/models"

os.chdir(log_file_path)
day = datetime.now().strftime('%d_%m_%Y')
logging.basicConfig(filename='extend_universal_for_gapfilling.log'.format(day), level=logging.INFO, filemode='w')
logger = logging.getLogger(__name__)

# universal reaction bag for model generation
os.chdir(model_path)
universal_model = cobra.io.read_sbml_model('universal_model_updated.xml')

if 'bof_c' in [r.id for r in universal_model.reactions]:
    universal_model.remove_reactions([universal_model.reactions.bof_c])
    universal_model.repair()

# extend universal by curated model
pf_model = cobra.io.read_sbml_model('iPfal19.xml')
len_univ_rxns = len(universal_model.reactions)
for rxn in pf_model.reactions:
    if rxn.id not in [r.id for r in universal_model.reactions]:
        if len(set(['hb_c','hb_e','hemozoin_c','hemozoin_e','hemozoin_fv']).intersection(set([met.id for met in rxn.metabolites.keys()]))) == 0:
            mets = [x.metabolites for x in [rxn]]
            all_keys = set().union(*(d.keys() for d in mets))
            for key in all_keys:
                if key.id not in [m.id for m in universal_model.metabolites]:
                    universal_model.add_metabolites([key.copy()])
            universal_model.add_reactions([rxn.copy()]) # extend universal by Pf reactions, but remove gene IDs
            if len(rxn.genes) > 0:
                genes = rxn.genes
                universal_model.reactions.get_by_id(rxn.id).gene_reaction_rule = ''
                for gene in genes:
                    if len(universal_model.genes.get_by_id(gene.id).reactions) == 0:
                           gene_list = [universal_model.genes.get_by_id(gene.id)]
                           cobra.manipulation.delete.remove_genes(universal_model, gene_list, remove_reactions=False)
    else: #rxn.id IS in [r.id for r in universal_model but in PF its reversible and in universal its irreversible, make reversible
        if rxn.lower_bound < universal_model.reactions.get_by_id(rxn.id).lower_bound:
            universal_model.reactions.get_by_id(rxn.id).lower_bound = rxn.lower_bound
        if rxn.upper_bound > universal_model.reactions.get_by_id(rxn.id).upper_bound:
            universal_model.reactions.get_by_id(rxn.id).upper_bound = rxn.upper_bound
if len(universal_model.reactions) <= len_univ_rxns:
    logger.info('ERROR - universal model does not have Pf reactions added!')
universal_model.repair()


#for rxn in pf_model.reactions:
#    if rxn.id in [r.id for r in universal_model.reactions] and rxn.reaction != universal_model.reactions.get_by_id(rxn.id).reaction:
#        logger.info('{} is in both universal and iPfal19 models, but with a different reaction string.'.format(rxn.id))
#        logger.info('universal:')
#        logger.info(universal_model.reactions.get_by_id(rxn.id).reaction)
#        logger.info('iPfal19:')
#        logger.info(rxn.reaction)

# extend by reactions in all models
# this is essential because we moved reactions in the universal model into other compartments when building each model
# and then (for gapfilling) prune the universal model to contain only the reaction relevant to that species (ie only the
# reactions that are in the right compartments)
os.chdir(model_path)
for file in glob.glob("final_denovo_*.json"):
    add_model = cobra.io.load_json_model(file)
    for rxn in add_model.reactions:
        if rxn.id not in [r.id for r in universal_model.reactions]:
            mets = [x.metabolites for x in [rxn]]            
            all_keys = set().union(*(d.keys() for d in mets))
            for key in all_keys:
                if key.id not in [m.id for m in universal_model.metabolites]:
       	            universal_model.add_metabolites([key.copy()])
            universal_model.add_reactions([rxn.copy()])
            if not (rxn.gene_reaction_rule == ''):     
                universal_model.reactions.get_by_id(rxn.id).gene_reaction_rule = ''
    cobra.manipulation.delete.remove_genes(universal_model, universal_model.genes, remove_reactions=False)
    universal_model.repair()
    logger.info('extended by_{}'.format(file))

# remove Biomass reactions
rxn_list_to_delete = [r.id for r in universal_model.reactions if r.id.startswith('BIOMASS_')]
universal_model.remove_reactions(rxn_list_to_delete)
rxn_list_to_delete = [r.id for r in universal_model.reactions if r.id.startswith('biomass')]
universal_model.remove_reactions(rxn_list_to_delete)
rxn_list_to_delete = [r.id for r in universal_model.reactions if r.id.startswith('generic_biomass')]
universal_model.remove_reactions(rxn_list_to_delete)
universal_model.repair()
logger.info('removed biomasses from universal')

# add exchange reactions
for met in universal_model.metabolites:
    if met.id.endswith('_e'):
        if 'EX_'+met.id not in universal_model.reactions:
            universal_model.add_boundary(met,type = "exchange")

# tidy
universal_model = hf.update_universal_model(universal_model)
universal_model, unused = hf.prune_unused_metabolites2(universal_model)
universal_model.repair()

# make sure all reactions can carry flux, problem with some versions of the universal model
for rxn in universal_model.reactions:
    if rxn.lower_bound == 0 and rxn.upper_bound == 0:
        logging.info(rxn.id + ' has bounds == 0 in '+key)
        rxn.lower_bound = -1000.
        rxn.upper_bound = 1000.
        # NOTHING SHOULD PRINT - this was a problem in CarveMe

os.chdir(model_path)
cobra.io.save_json_model(universal_model, "extended_universal_model_for_gapfilling.json")
cobra.io.write_sbml_model(universal_model, "extended_universal_model_for_gapfilling.xml")
