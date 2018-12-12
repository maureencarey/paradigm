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
from datetime import datetime

data_path = "/scratch/mac9jc/paradigm/data"
model_path = "/scratch/mac9jc/paradigm/models"
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
if 'with_biomass_denovo_' in SPECIES_ID:
    SPECIES_ID = SPECIES_ID.split('with_biomass_denovo_')[1]
if 'ortho_' in SPECIES_ID:
    SPECIES_ID = SPECIES_ID.split('ortho_')[1]
if 'DIY' in SPECIES_ID:
    SPECIES_ID = SPECIES_ID[5:]
if SPECIES_ID == 'iPfal18':
    SPECIES_ID = 'Pfalciparum3D7'

day = datetime.now().strftime('%d_%m_%Y')
logging.basicConfig(filename='step6_{}_{}.log'.format(SPECIES_ID,day), level=logging.INFO, filemode='w')
logger = logging.getLogger(__name__)
logger.info('BEGIN STEP 6')
    
# modified for Rivanna: read in the models
model_dict = {}
model_dict[SPECIES_ID] = cobra.io.load_json_model(model_fname)
logger.info('loaded model')

for species, model in model_dict.items():
    if 'biomass' not in [r.id for r in model.reactions]:
        logger.info('biomass not in reactions anymore')

os.chdir(model_path)
universal_model = cobra.io.load_json_model('universal_model_oct26_2018.json')

# extend universal by curated model
pf_model = cobra.io.load_json_model('iPfal18.json')
for rxn in pf_model.reactions:
    if rxn.id not in [r.id for r in universal_model.reactions]:
        if len(set(['hb_c','hb_e','hemozoin_c','hemozoin_e','hemozoin_fv']).intersection(set([met.id for met in rxn.metabolites.keys()]))) == 0:
            universal_model.add_reactions([rxn.copy()])

logger.info('loaded universal')
genome_ids = pd.read_csv("auxotrophies_mapping_to_genomeID.csv",header = None).T
new_header = genome_ids.iloc[0]
genome_ids = genome_ids[1:]
genome_ids.columns = new_header

met_ids = pd.read_csv("auxotrophies_mapping_to_metID.csv")

os.chdir(data_path)
gapfilling_tasks = pd.read_excel("auxotrophies_references.xlsx",skiprows = 1)
idx2 = met_ids['BiGG'].notnull()
temp_met_ids = met_ids.loc[idx2]
met_dict = pd.Series(temp_met_ids.BiGG.values, index = temp_met_ids.Metabolite).to_dict()

# change met ids to model met ids if possible
for x in gapfilling_tasks.index:
    if gapfilling_tasks['Metabolite'].iloc[x] in met_dict.keys():
        gapfilling_tasks['Metabolite'].iloc[x] = met_dict[gapfilling_tasks['Metabolite'].\
        iloc[x]]

for species, model in model_dict.items():
    logger.info(species)
    logger.info(model.objective.expression)
    logger.info('')

tasks_dict = dict()
for species, model in model_dict.items():
    if sum(genome_ids['strain_ID'].isin([species])): # does the model have any tasks?
        idx = genome_ids.index[genome_ids['strain_ID'] == species].tolist()
        sp = genome_ids.loc[idx]['species'] # get species name
        species_string = sp[sp.index[0]]
        
        # tasks for that species
        gf_species = gapfilling_tasks[["Metabolite",species_string]]
        indx = gf_species[species_string].notnull()
        if sum(indx) >0:
            gf_species.loc[indx]
            
            # types of tasks
            consumption_species = gf_species[gf_species[species_string].str.contains(\
            "uptake|essential|both|rescue")==True]
            production_species =gf_species[gf_species[species_string].str.contains(\
            "prod|both")==True]
            
            tasks_dict[species] = {"consumption":consumption_species['Metabolite'],
                "production":production_species['Metabolite']}
    if 'biomass' not in [r.id for r in model.reactions]:
        logger.info('biomass not in reactions anymore')
    if 'hb_c' in [m.id for m in model.metabolites]:
        logger.info('HEMOGLOBIN PRESENT')

# tasks_dict # key = species, value = dictionary
# tasks_dict[species][consumption] = dataframe of met Ids with column name 'Metabolite'
# tasks_dict[species][production] = dataframe of met Ids with column name 'Metabolite'

from cobra.flux_analysis.parsimonious import add_pfba

def flatten_mixed_list(list_of_interest):
    new_list = list()
    for x in list_of_interest:
        if isinstance(x,list):
            new_list.extend(x)
        else:
            new_list.append(x)
    return(new_list)

# pFBA based gapfilling, implementation from Greg Medlock
def pfba_gapfill_implementation(input_model, universal_model_ex, objective_reaction_id):
    # objective_reaction is a reaction id
    
    universal_model_pfba = universal_model_ex.copy()
    add_pfba(universal_model_pfba)
    
    # penalize adding demand reactions
    coef = universal_model_pfba.objective.get_linear_coefficients(universal_model_pfba.variables)
    for key,value in coef.items():
        if key.name.startswith('DM_') or key.name.startswith('SK_'):
            coef[key] = 1000.
        elif key.name.startswith('EX_'):
            coef[key] = 50.
    for x in universal_model_pfba.variables:
        if x.name.startswith('DM_') or x.name.startswith('SK_') or x.name.startswith('EX_'):
            universal_model_pfba.objective.set_linear_coefficients(coef)

    rxns_to_remove = [rxn for rxn in universal_model_pfba.reactions if rxn.id \
                      in [rxn.id for rxn in input_model.reactions]]
    universal_model_pfba.remove_reactions(rxns_to_remove)
    universal_model_pfba.add_reactions([rxn for rxn in input_model.reactions])
    universal_model_pfba.reactions.get_by_id(objective_reaction_id).lower_bound = 0.05

    solution = universal_model_pfba.optimize()
    if solution.status == 'infeasible':
        logging.info('pFBA gapfilling for {} is infeasible!'.format(objective_reaction_id))
    else:
        logging.info('pFBA gapfilling - feasible')
    get_fluxes = set([r.id for r in universal_model_pfba.reactions]) - set([rxn.id for rxn in input_model.reactions]) ### ERROR IN ORIGINAL CODE HAD MODEL instead of INPUT_MODEL
    add_reactions_to_model = [rxn for rxn in get_fluxes if abs(solution.x_dict[rxn]) > 1E-8]

    # double check
    logging.info('double checking pFBA solution')
    test_model = input_model.copy()
    add_reactions_list = [universal_model_pfba.reactions.get_by_id(r).copy() for r in add_reactions_to_model]
    test_model.add_reactions(add_reactions_list)
    sol = test_model.optimize()
    if sol.status == 'infeasible':
        logging.info('pFBA solution is infeasible!')
    return(add_reactions_to_model)

for species, model in model_dict.items():
    
    logger.info(species)
    model.solver = 'gurobi'
    if 'biomass' in [r.id for r in model.reactions]:
        model.objective = 'biomass'
        logger.info('was able to set objective')
    add_reactions_list = list()
    add_mets_list = list()
    
    mets_to_prod = tasks_dict[species]["production"]
    mets_to_consume = tasks_dict[species]["consumption"]
    produce_met = list()
    consume_met = list()
    all_mets = list()
    cannot_gapfill = list()
    for met in mets_to_prod:
        if met+'_c' in universal_model.metabolites:
            produce_met.append(met)
            all_mets.append(met)
        else:
            cannot_gapfill.append(met)
    for met in mets_to_consume:
        if met+'_c' in universal_model.metabolites:
            consume_met.append(met)
            all_mets.append(met)
        else:
            cannot_gapfill.append(met)
    all_mets = list(set(all_mets))

    if 'hb_c' in [m.id for m in model.metabolites]:
        logger.info('HEMOGLOBIN PRESENT2')
    for met in all_mets:
        logger.info('')
        logger.info(met)
        extracel_met = met+'_e'
        cyto_met = met+'_c'
    
        gf_model = model.copy()
        gf_universal = universal_model.copy()
        # make sure you are making any mets by reversing biomass production
        if 'biomass' in [r.id for r in gf_model.reactions]:
            gf_model.reactions.get_by_id('biomass').lb = 0.
            gf_model.reactions.get_by_id('biomass').ub = 0.
        elif 'generic_biomass' in [r.id for r in gf_model.reactions]:
            gf_model.reactions.get_by_id('generic_biomass').lb = 0.
            gf_model.reactions.get_by_id('generic_biomass').ub = 0.
    
        # add metabolites if necessary
        if cyto_met not in gf_model.metabolites:
            # already filtered list for only mets_c in universal
            logger.info('add cyto met')
            gf_model.add_metabolites([gf_universal.metabolites.get_by_id(cyto_met).\
            copy()])
            add_mets_list.append(gf_universal.metabolites.get_by_id(cyto_met).\
            copy())
        if extracel_met not in gf_model.metabolites:
            logger.info('add extracel met')
            if extracel_met in gf_universal.metabolites:
                gf_model.add_metabolites([gf_universal.metabolites.get_by_id(\
                extracel_met).copy()])
                add_mets_list.append(gf_universal.metabolites.get_by_id(\
                extracel_met).copy())
            else:
                met = Metabolite(extracel_met,
                                formula = gf_universal.metabolites.get_by_id(\
                                cyto_met).formula,
                                name = gf_universal.metabolites.get_by_id(\
                                cyto_met).name,
                                compartment='e')
                gf_model.add_metabolites([met])
                add_mets_list.append(met)

        # add objective reaction
        gf_model.add_boundary(gf_model.metabolites.get_by_id(cyto_met), type = "demand")
        gf_model.objective = 'DM_'+cyto_met

        # add exchange
        if 'EX_'+extracel_met not in gf_model.reactions:
            logger.info('add exchange')
            gf_model.add_boundary(gf_model.metabolites.get_by_id(extracel_met), \
            type="exchange")
            add_reactions_list.append(gf_model.reactions.get_by_id('EX_'+\
            extracel_met).copy())
        
        if met in list(produce_met):
            logger.info('produce')
            gf_model.reactions.get_by_id('EX_'+extracel_met).lower_bound = 0.
            gf_model.reactions.get_by_id('EX_'+extracel_met).upper_bound = 0.
            if model.slim_optimize() < 0.01: # gapfill if can't produce
                solution = pfba_gapfill_implementation(gf_model, gf_universal, \
                'DM_'+cyto_met)
                logger.info('solution')
                logger.info(solution)
                if len(solution) > 1:
                    for rxn_id in solution:
                        if rxn_id != 'DM_'+cyto_met:
                            add_reactions_list.append(gf_universal.reactions.\
                            get_by_id(rxn_id).copy())
                        else:
                            if solution[0] != 'DM_'+cyto_met:
                               add_reactions_list.append(gf_universal.reactions.\
                               get_by_id(solution[0]).copy())
            else:
                logger.info('no gapfilling needed')

        if met in list(consume_met):
            logger.info('consume')
            gf_model.reactions.get_by_id('EX_'+extracel_met).lb = -1000.
            gf_model.reactions.get_by_id('EX_'+extracel_met).ub = 1000.
            t = False
            if gf_model.slim_optimize() < 0.01: # optimize and no growth
                # add import reaction (don't use gapfill because we want an import,
                # not biosynthetic solution)
                if extracel_met in gf_universal.metabolites:
                    extracellular_reactions = gf_universal.metabolites.get_by_id(\
                    extracel_met).reactions
                    intracellular_reactions = gf_universal.metabolites.get_by_id(\
                    cyto_met).reactions
                    possible_rxns = list(extracellular_reactions.intersection(\
                    intracellular_reactions))
                    if len(possible_rxns) >0:
                        if len(possible_rxns) > 1:
                            rxn_keep = ''
                            for rxn in possible_rxns:
                                if len(rxn.metabolites) == 2:
                                    # add reaction that has no other substrates
                                    rxn_keep = rxn
                                else:
                                    pass
                            if rxn_keep != '':
                                add_reactions_list.append(rxn_keep)
                            else: add_reactions_list.append(possible_rxns[0])
                            logger.info('possible_rxns')
                            logger.info(rxn_keep)
                        else:
                            add_reactions_list.append(possible_rxns[0])
                            logger.info('possible_rxns')
                            logger.info(possible_rxns[0])
                    
                    else:
                        t = True # if no transport in universal, make transport rxn
                else:
                    t = True # if extracel met isnt in universal, make transport rxn
            else:
                logger.info('no gapfilling needed')
                    # WORKS
         
            if t == True:
                reaction = Reaction(met.upper()+'t')
                reaction.name = met + 'transport'
                reaction.subsystem = 'Transport'
                reaction.lower_bound = -1000.
                reaction.upper_bound = 1000.
                reaction.add_metabolites({gf_model.metabolites.get_by_id(extracel_met): \
                -1.0, gf_model.metabolites.get_by_id(cyto_met): 1.0 })
                add_reactions_list.append(reaction)
                    
    if 'biomass' not in [r.id for r in gf_model.reactions]:
        logger.info('biomass not in reactions anymore')
    if 'hb_c' in [m.id for m in model.metabolites]:
        logger.info('HEMOGLOBIN PRESENT3')
# is it used in any reaction intracellularly? 
# if not, cannot do anything about it, except curate in future

    add_reactions_list = flatten_mixed_list(add_reactions_list)
    df = pd.DataFrame({'rxns_added':add_reactions_list})
    df.to_csv('gapfilling_additions_{0}_tasks.csv'.format(species))

    for met in add_mets_list:
        if met.id not in [m.id for m in model.metabolites]:
            model.add_metabolites([met])
    for rxn in add_reactions_list:
        if rxn.id not in [r.id for r in model.reactions]:
            model.add_reactions([rxn])
    if 'hb_c' in [m.id for m in model.metabolites]:
        logger.info('HEMOGLOBIN PRESENT4')

    logger.info("---------------------------------------------------")
    logger.info(species)
    logger.info('no. reactions added:')
    logger.info(len(set(add_reactions_list)))
    logger.info('no. metabolites added:')
    logger.info(len(set(add_mets_list)))
    # ADD FIRST
    model_dict[species] = model
    # save models

GENERIC_BIOMASS = 'generic_biomass'
P_BIOMASS = 'biomass'
for species, model in model_dict.items():

    logger.info(species)
    model.solver = 'gurobi'
    gf_mod_list1 = list()
    gf_mod_list2 = list()
    if 'biomass' in [r.id for r in model.reactions]:
        model.reactions.get_by_id('biomass').lb = 0.
        model.reactions.get_by_id('biomass').ub = 1000.
    elif 'generic_biomass' in [r.id for r in model.reactions]:
        model.reactions.get_by_id('generic_biomass').lb = 0.
        model.reactions.get_by_id('generic_biomass').ub = 1000.

    if GENERIC_BIOMASS in [r.id for r in model.reactions]:
        gf_model = model.copy()
        gf_universal = universal_model.copy()
        gf_model.objective = GENERIC_BIOMASS
        if gf_model.slim_optimize() < 0.01: # gapfill if can't produce
             solution = pfba_gapfill_implementation(gf_model, gf_universal, \
             GENERIC_BIOMASS)
             if len(solution) >= 1:
                 for rxn_id in solution:
                     add_reactions_list.append(gf_universal.reactions.\
                     get_by_id(rxn_id).copy())
                     gf_mod_list1.append(gf_universal.reactions.\
                     get_by_id(rxn_id).copy())
        os.chdir(data_path)
        df = pd.DataFrame({'rxns_added':gf_mod_list1})
        df.to_csv('gapfilling_additions_{0}_generic_biomass.csv'.format(SPECIES_ID_old))
        logger.info("wrote generic bioamss file")
    else:
        logger.info('error: no generic biomass reaction')

    if 'hb_c' in [m.id for m in model.metabolites]:
        logger.info('HEMOGLOBIN PRESENT5')
    
    logger.info('biomass in reactions')
    logger.info('biomass' in [r.id for r in model.reactions])
#    logger.info(model.objective.expression)

    if 'biomass' in [r.id for r in model.reactions]:
        gf_model = model.copy()
        gf_universal = universal_model.copy()
        gf_model.objective = 'biomass'
        logger.info('set objective')
        logger.info(gf_model.slim_optimize())
        if gf_model.slim_optimize() < 0.01: # gapfill if can't produce
             solution = pfba_gapfill_implementation(gf_model, gf_universal, \
             'biomass')
             if len(solution) >= 1:
                 logger.info('reactions to add')
                 logger.info(len(solution))
                 for rxn_id in solution:
                     add_reactions_list.append(gf_universal.reactions.\
                     get_by_id(rxn_id).copy())
                     gf_mod_list2.append(gf_universal.reactions.\
                     get_by_id(rxn_id).copy())
        os.chdir(data_path)
        df = pd.DataFrame({'rxns_added':gf_mod_list2})
        df.to_csv('gapfilling_additions_{0}_species_biomass.csv'.format(SPECIES_ID_old))
        logger.info("wrote species biomass file")

    if 'hb_c' in [m.id for m in model.metabolites]:
        logger.info('HEMOGLOBIN PRESENT6')

    add_reactions_list2 = list(set(flatten_mixed_list([gf_mod_list1,gf_mod_list2])))
    logger.info('going to add these reactions:')
    logger.info(add_reactions_list2)
    for rxn in add_reactions_list2:
        if rxn.id not in [r.id for r in model.reactions]:
            model.add_reactions([rxn])
    logger.info('added reactions')

    if 'hb_c' in [m.id for m in model.metabolites]:
        logger.info('HEMOGLOBIN PRESENT7')
    cobra.io.save_json_model(model, 'gf_'+SPECIES_ID_old+'.json')
    
    
