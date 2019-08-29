# Input with_biomass_speciesname.json and iPfal19.json

import cobra
import pandas as pd
import os
import subprocess
import glob
import json
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis.parsimonious import add_pfba
#import helper_functions_2 as hf2
import sys
sys.path.append(os.path.abspath("/home/mac9jc/paradigm/"))
import helper_functions as hf
import argparse
import logging
from datetime import datetime
from itertools import chain
from cobra.util.solver import linear_reaction_coefficients
from optlang.interface import OPTIMAL

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
SPECIES_ID_old = SPECIES_ID
if 'denovo_' in SPECIES_ID:
    SPECIES_ID = SPECIES_ID.split('denovo_')[1]
if 'with_biomass_' in SPECIES_ID:
    SPECIES_ID = SPECIES_ID.split('with_biomass_')[1]
if 'ortho_' in SPECIES_ID:
    SPECIES_ID = SPECIES_ID.split('ortho_')[1]
if 'DIY' in SPECIES_ID:
    SPECIES_ID = SPECIES_ID[5:]
if SPECIES_ID == 'iPfal19':
    SPECIES_ID = 'Pfalciparum3D7'

day = datetime.now().strftime('%d_%m_%Y')
logging.basicConfig(filename='step6_no_ortho_{}_{}.log'.format(SPECIES_ID_old,day), level=logging.INFO, filemode='w')
logger = logging.getLogger(__name__)
logger.info('BEGIN STEP 6 for a Plasmodium model - WITHOUT ORTHOLOGY TRANSFORMATION')

os.chdir(model_path)
model = cobra.io.load_json_model(model_fname)
model.repair()
model.solver = 'glpk'
logger.info('loaded model')

def validate(original_model, reactions):
    with original_model as model:
        mets = [x.metabolites for x in reactions]
        all_keys = set().union(*(d.keys() for d in mets))
        for key in all_keys:
            if key.id not in [m.id for m in model.metabolites]:
                model.add_metabolites([key.copy()])
        model.add_reactions(reactions)
        model.slim_optimize()
        return (model.solver.status == OPTIMAL and model.solver.objective.value >= 0.0001)

# write here because it uses logger
# pFBA based gapfilling, implementation from Greg Medlock in his Medusa package
def pfba_gapfill_implementation(input_model, universal_model_ex, objective_reaction_id):
    # objective_reaction is a reaction id

    gapfiller = universal_model_ex.copy()
    if linear_reaction_coefficients(gapfiller):
        logger.info('ERROR: universal model has an existing objective rxn')

    # get the original objective from the model being gapfilled
    model_to_gapfill = input_model.copy()
    original_objective = linear_reaction_coefficients(model_to_gapfill)
    # convert to IDs to avoid issues with model membership when these reactions
    # are added to gapfiller
    original_objective = {rxn.id:original_objective[rxn] for rxn
        in original_objective.keys()}

    # get the reactions in the original model, which need to be removed from
    # the universal if present. This cannot catch identical reactions that do
    # not share IDs, so make sure your model and universal are in the same
    # namespace.
    rxns_to_remove = [rxn for rxn in gapfiller.reactions if rxn.id in \
                      [rxn.id for rxn in model_to_gapfill.reactions]]
    gapfiller.remove_reactions(rxns_to_remove)
    gapfiller.repair()

    # get the list of reactions currently in the gapfiller, which are the ones
    # we will need to check for flux after solving the problem (e.g. these are
    # the reactions we are considering adding to the model)
    get_fluxes = [rxn.id for rxn in gapfiller.reactions]

    # add the reactions from the model to the gapfiller, which are not
    # included in the pFBA formulation, and thus flux is not penalized
    # through them.
    original_model_reactions = [rxn.copy() for rxn in model_to_gapfill.reactions]
    mets = [x.metabolites for x in original_model_reactions]
    all_keys = set().union(*(d.keys() for d in mets))
    for key in all_keys:
        if key.id not in [m.id for m in gapfiller.metabolites]:
            gapfiller.add_metabolites([key.copy()])
    gapfiller.add_reactions(original_model_reactions)
    original_reaction_ids = [reaction.id for reaction in original_model_reactions]

    # Add the pFBA constraints and objective (minimizes sum of fluxes)
    add_pfba(gapfiller)

    # ORDER MATTERS HERE: penalize adding exchange reactions
    coefficients = (gapfiller.objective.get_linear_coefficients(gapfiller.variables))
    penalize_these_EX_reactions = [r.id for r in gapfiller.reactions if r.id.startswith('EX_')]
    reaction_variables_ex = (((gapfiller.reactions.get_by_id(reaction).forward_variable),
                              (gapfiller.reactions.get_by_id(reaction).reverse_variable))
                             for reaction in penalize_these_EX_reactions)
    variables_ex = chain(*reaction_variables_ex)
    for variable in variables_ex: coefficients[variable] = 10.0

    # penalize adding pseudoreactions like demand and sink reactions
    penalize_these_PS_reactions = [r.id for r in gapfiller.reactions\
                                   if r.id.startswith('DM_') or r.id.startswith('SK_')]
    reaction_variables_ps = (((gapfiller.reactions.get_by_id(reaction).forward_variable),
                              (gapfiller.reactions.get_by_id(reaction).reverse_variable))
                             for reaction in penalize_these_PS_reactions)
    variables_ps = chain(*reaction_variables_ps)
    for variable in variables_ps: coefficients[variable] = 100.0

    # set the linear coefficients for reactions in the original model to 0
    reaction_variables = (((gapfiller.reactions.get_by_id(reaction).forward_variable),
                           (gapfiller.reactions.get_by_id(reaction).reverse_variable))
                           for reaction in original_reaction_ids)
    variables = chain(*reaction_variables)
    for variable in variables: coefficients[variable] = 0.0
    gapfiller.objective.set_linear_coefficients(coefficients)

    # set a constraint on flux through the original objective
    for reaction in original_objective.keys():
        #logger.info({'this is the objective, is it what i expect?':gapfiller.reactions.get_by_id(reaction).id})
        gapfiller.reactions.get_by_id(reaction).lower_bound = 0.01
    
    # get solution
    solution = gapfiller.optimize()
    filtered_solution = {rxn:solution.fluxes[rxn] for rxn in\
        get_fluxes if abs(solution.fluxes[rxn]) > 1E-10}
    add_rxns = [universal_model_ex.reactions.get_by_id(rxn).copy() for \
                rxn in filtered_solution.keys()]
    cycle_reactions = set([rxn.id for rxn in add_rxns])

    # validate that the proposed solution restores flux through the
    # objective in the original model
    # set the bounds on the original model to represent media
    # and validate the gapfill solution
    if not validate(model_to_gapfill, [universal_model_ex.reactions.get_by_id(rxn_id).copy() for rxn_id in cycle_reactions]):
        logger.info('INFEASBIBLE gapfill solution for '+objective_reaction_id)
        #raise RuntimeError('Failed to validate gapfilled model, '
        #                              'try lowering the flux_cutoff through '
        #                               'inclusion_threshold')
    return cycle_reactions

os.chdir(model_path)
universal_model = cobra.io.read_sbml_model("extended_universal_model_for_gapfilling.xml")
iPfal19 = cobra.io.read_sbml_model("iPfal19.xml")
universal_model.repair()
universal_model.solver = 'glpk'
iPfal19.repair()
iPfal19.solver = 'glpk'

## prep for removing universal reactions that are in the wrong compartment for this model
# database mapping - this is for release 42, must update if using a different EuPathDB release
plasmodb = ["PadleriG01","PbergheiANKA","PbillcollinsiG01","PblacklockiG01","Pchabaudichabaudi","PcoatneyiHackeri","PcynomolgiB","PcynomolgiM","Pfalciparum3D7","Pfalciparum7G8","PfalciparumCD01","PfalciparumDd2","PfalciparumGA01","PfalciparumGB4","PfalciparumGN01","PfalciparumHB3","PfalciparumIT","PfalciparumKE01","PfalciparumKH01","PfalciparumKH02","PfalciparumML01","PfalciparumSD01","PfalciparumSN01","PfalciparumTG01","PfragileNilgiri","PgaboniG01","PgaboniSY75","Pgallinaceum8A","PinuiSanAntonio1","PknowlesiH","PknowlesiMalayanPk1A","PmalariaeUG01","PovalecurtisiGH01","PpraefalciparumG01","PreichenowiCDC","PreichenowiG01","PrelictumSGS1-like","PvinckeipetteriCR","Pvinckeivinckeivinckei","Pvivax-likePvl01","PvivaxP01","PvivaxSal1","Pyoeliiyoelii17X","Pyoeliiyoelii17XNL","PyoeliiyoeliiYM"]

if SPECIES_ID in plasmodb: # Plasmodium = cytosol, extracellular, mitochondrdia, apicoplast, food vacuole
    compartments = ["_c","_e","_m","_ap","_fv"]
else:
    logger.info('NOT A PLASMODIUM')
    compartments = ["_c","_e"]

logger.info('reading experimental data')
# experimental data
os.chdir(data_path)
genome_ids = pd.read_csv("auxotrophies_mapping_to_genomeID.csv",header = None).T
new_header = genome_ids.iloc[0]
genome_ids = genome_ids[1:]
genome_ids.columns = new_header
met_ids = pd.read_csv("auxotrophies_mapping_to_metID.csv")
gapfilling_tasks = pd.read_excel("auxotrophies_references.xlsx",skiprows = 1)
idx2 = met_ids['BiGG'].notnull()
temp_met_ids = met_ids.loc[met_ids['BiGG'].notnull()].copy()
met_dict = pd.Series(met_ids.loc[idx2].BiGG.values, index = met_ids.loc[idx2].Metabolite).to_dict()

# change met ids to model met ids if possible
for x in gapfilling_tasks.index:
    if gapfilling_tasks['Metabolite'].iloc[x] in met_dict.keys():
        gapfilling_tasks['Metabolite'].iloc[x] = met_dict[gapfilling_tasks['Metabolite'].iloc[x]]

if SPECIES_ID in plasmodb:
    r_len = len(universal_model.reactions)
    for rxn_id in ['EX_hb_e','HBtr','HMGLB','HMBZ','HMBZex','EX_hemozoin_e']:
        rxn = iPfal19.reactions.get_by_id(rxn_id)
        mets = [x.metabolites for x in [rxn]]
        all_keys = set().union(*(d.keys() for d in mets))
        for key in all_keys:
            if key.id not in [m.id for m in universal_model.metabolites]:
                universal_model.add_metabolites([key.copy()])
        universal_model.add_reactions([rxn.copy()])
        if len(rxn.genes) > 0:
            genes = rxn.genes
            universal_model.reactions.get_by_id(rxn.id).gene_reaction_rule = ''
            for gene in genes:
                if len(universal_model.genes.get_by_id(gene.id).reactions) == 0:
                    gene_list = [universal_model.genes.get_by_id(gene.id)]
                    cobra.manipulation.delete.remove_genes(universal_model, gene_list, remove_reactions=False)
    if len(universal_model.reactions) <= r_len:
        logger.info('error in adding reactions and removing genes')
    for rxn_id in ['EX_hb_e','HBtr','HMGLB','HMBZ','HMBZex','EX_hemozoin_e']:
        if rxn_id not in [r.id for r in universal_model.reactions]:
            logger.info('{} not correctly added to universal for gapfilling'.format(rxn_id))
else: logger.info('not in plasmodb, thus HB reactions were not added')
universal_model.repair()

mets_to_prod = list()
mets_to_consume = list()
if sum(genome_ids['strain_ID'].isin([SPECIES_ID])): # does the model have any tasks?
    idx = genome_ids.index[genome_ids['strain_ID'] == SPECIES_ID].tolist()
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
        
        mets_to_prod = production_species['Metabolite'].tolist()
        mets_to_consume = consumption_species['Metabolite'].tolist()


if 'biomass' in [r.id for r in model.reactions]:
    model.objective = 'biomass'
    logger.info('was able to set objective')
add_reactions_list = list()
add_mets_list = list()

# remove reactions that are not in eligible compartments
universal_model_for_species = universal_model.copy()
for rxn in universal_model_for_species.reactions:
    rxn_metabolite_list = [met.id for met in rxn.metabolites]
    # CHECK set(metabolite.compartment) # check if attribute is not blank
    rxn_metabolite_comp_list = [hf.get_comp(universal_model_for_species,m) for m in rxn_metabolite_list]
    if len(hf.unaccept_comp_intersection(rxn_metabolite_comp_list, compartments))>0:
        universal_model_for_species.remove_reactions([rxn])
universal_model_for_species, unused  = cobra.manipulation.delete.prune_unused_metabolites(universal_model_for_species)
universal_model_for_species.repair()

produce_met = list()
consume_met = list()
all_mets = list()
cannot_gapfill = list()
if mets_to_prod:
    for met in mets_to_prod:
        if met+'_c' in universal_model.metabolites:
            produce_met.append(met)
            all_mets.append(met)
        else:
            cannot_gapfill.append(met)
if mets_to_consume:
    for met in mets_to_consume:
        if met+'_c' in universal_model.metabolites:
            consume_met.append(met)
            all_mets.append(met)
        else:
            cannot_gapfill.append(met)
all_mets = list(set(all_mets))

for met in all_mets:
    logger.info('')
    logger.info(met)
    extracel_met = met+'_e'
    cyto_met = met+'_c'

    gf_model = model.copy()
    if 'biomass' in [r.id for r in gf_model.reactions]:
       	gf_model.remove_reactions([gf_model.reactions.get_by_id('biomass')])
    if 'generic_biomass' in [r.id for r in gf_model.reactions]:
        gf_model.remove_reactions([gf_model.reactions.get_by_id('generic_biomass')])
    gf_universal = universal_model_for_species.copy()

    # add metabolites if necessary
    if cyto_met not in gf_model.metabolites:
        # already filtered list for only mets_c in universal
        logger.info('add cyto met')
        gf_model.add_metabolites([gf_universal.metabolites.get_by_id(cyto_met).copy()])
        add_mets_list.append(gf_universal.metabolites.get_by_id(cyto_met).copy())
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
    if 'DM_'+cyto_met not in [r.id for r in gf_model.reactions]:
        gf_model.add_boundary(gf_model.metabolites.get_by_id(cyto_met), type = "demand")
    gf_model.objective = 'DM_'+cyto_met

    # add exchange
    if 'EX_'+extracel_met not in gf_model.reactions:
        gf_model.add_boundary(gf_model.metabolites.get_by_id(extracel_met), type="exchange")
        add_reactions_list.append(gf_model.reactions.get_by_id('EX_'+extracel_met).copy())

    if met in list(produce_met):
        gf_model.reactions.get_by_id('EX_'+extracel_met).lower_bound = 0.
        gf_model.reactions.get_by_id('EX_'+extracel_met).upper_bound = 0.
        if not (gf_model.slim_optimize() > 0.01): # gapfill if can't produce
            # is gapfilling even possible?
            with gf_universal as model:
                for rxn in gf_model.reactions:
                    if rxn.id in [r.id for r in model.reactions]:
                        model.remove_reactions([model.reactions.get_by_id(rxn.id)])
                    mets = [x.metabolites for x in [rxn]]
                    all_keys = set().union(*(d.keys() for d in mets))
                    for key in all_keys:
                        if key.id not in [m.id for m in model.metabolites]:
                            model.add_metabolites([key.copy()])
                    model.add_reactions([rxn.copy()])
                model.repair()
                model.objective = 'DM_'+cyto_met
                if model.slim_optimize() < 0.01:
                    logger.info('INFEASBILE: gapfilling is not possible for production of '+extracel_met)
        
            # if yes, gapfill
            solution = pfba_gapfill_implementation(gf_model, gf_universal, 'DM_'+cyto_met)
            logger.info('solution to produce:')
            logger.info(solution)
            if len(solution) > 1:
                for rxn_id in solution:
                    add_reactions_list.append(gf_universal.reactions.get_by_id(rxn_id).copy())
        else:
            logger.info('no gapfilling needed')

    if met in list(consume_met):
        gf_model.reactions.get_by_id('EX_'+extracel_met).lb = -1000.
        gf_model.reactions.get_by_id('EX_'+extracel_met).ub = 1000.
        t = False
        if not (gf_model.slim_optimize() > 0.01): # optimize and no growth
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
                            if len(rxn.metabolites) == 2: rxn_keep = rxn # add if no other substrates
                            else: pass
                        if rxn_keep != '':
                            add_reactions_list.append(rxn_keep)
                        else: add_reactions_list.append(possible_rxns[0])
                        logger.info('solution to consume')
                        logger.info(rxn_keep)
                    else:
                        add_reactions_list.append(possible_rxns[0])
                        logger.info('solution to consume')
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
    
# is it used in any reaction intracellularly?
# if not, cannot do anything about it, except curate in future

add_reactions_list = hf.flatten_mixed_list(add_reactions_list)
df = pd.DataFrame({'rxns_added':add_reactions_list})
df.to_csv('gapfilling_additions_no_ortho_{0}_tasks.csv'.format(SPECIES_ID))

for met in add_mets_list: # these are copies, GOOD CHECK
    if met.id not in [m.id for m in model.metabolites]:
        model.add_metabolites([met])
for rxn in add_reactions_list:
    if rxn.id not in [r.id for r in model.reactions]:
        mets = [x.metabolites for x in [rxn]]
        all_keys = set().union(*(d.keys() for d in mets))
        for key in all_keys:
            if key.id not in [m.id for m in model.metabolites]:
                model.add_metabolites([key.copy()])
        model.add_reactions([rxn.copy()])

logger.info("---------------------------------------------------")
logger.info('no. reactions added:')
logger.info(len(set(add_reactions_list)))
logger.info('no. metabolites added:')
logger.info(len(set(add_mets_list)))

gf_mod_list1 = list()
gf_mod_list2 = list()
os.chdir(data_path)
if 'generic_biomass' in [r.id for r in model.reactions]:
    gf_model = model.copy()
    gf_model.reactions.get_by_id('generic_biomass').lb = 0.
    gf_model.reactions.get_by_id('generic_biomass').ub = 1000.
    if 'biomass' in [r.id for r in gf_model.reactions]:
       	gf_model.remove_reactions([gf_model.reactions.get_by_id('biomass')])
    gf_universal = universal_model_for_species.copy()

    gf_model.objective = 'generic_biomass'
    if not (gf_model.slim_optimize() > 0.01): # gapfill if can't produce
        # is gapfilling even possible?
        with gf_universal as model:
            for rxn in gf_model.reactions:
                if rxn.id in [r.id for r in model.reactions]:
                    model.remove_reactions([model.reactions.get_by_id(rxn.id)])
                mets = [x.metabolites for x in [rxn]]
                all_keys = set().union(*(d.keys() for d in mets))
                for key in all_keys:
                    if key.id not in [m.id for m in model.metabolites]:
                        model.add_metabolites([key.copy()])
                model.add_reactions([rxn.copy()])
            model.repair()
            model.objective = 'generic_biomass'
            if model.slim_optimize() < 0.01:
                logger.info('INFEASBILE: gapfilling is not possible for generic biomass')
    
        # if yes, gapfill
        solution = pfba_gapfill_implementation(gf_model, gf_universal, 'generic_biomass')
        if len(solution) >= 1:
            for rxn_id in solution:
                add_reactions_list.append(gf_universal.reactions.get_by_id(rxn_id).copy())
                gf_mod_list1.append(gf_universal.reactions.get_by_id(rxn_id).copy())
    df = pd.DataFrame({'rxns_added':gf_mod_list1})
    df.to_csv('gapfilling_additions_no_ortho_{0}_generic_biomass.csv'.format(SPECIES_ID_old))
    logger.info("wrote generic bioamss file")
else:
    logger.info('error: no generic biomass reaction')

if 'biomass' in [r.id for r in model.reactions]:
    gf_model = model.copy()
    gf_model.reactions.get_by_id('biomass').lb = 0.
    gf_model.reactions.get_by_id('biomass').ub = 1000.
    if 'generic_biomass' in [r.id for r in gf_model.reactions]:
        gf_model.remove_reactions([gf_model.reactions.get_by_id('generic_biomass')])
    gf_universal = universal_model_for_species.copy()

    gf_model.objective = 'biomass'
    if not (gf_model.slim_optimize() > 0.01): # gapfill if can't produce
        # is gapfilling even possible?
        with gf_universal as model:
            for rxn in gf_model.reactions:
                if rxn.id in [r.id for r in model.reactions]:
                    model.remove_reactions([model.reactions.get_by_id(rxn.id)])
                mets = [x.metabolites for x in [rxn]]
                all_keys = set().union(*(d.keys() for d in mets))
                for key in all_keys:
                    if key.id not in [m.id for m in model.metabolites]:
                        model.add_metabolites([key.copy()])
                model.add_reactions([rxn.copy()])
            model.repair()
            model.objective = 'biomass'
            if model.slim_optimize() < 0.01:
                logger.info('INFEASBILE: gapfilling is not possible for specific biomass')
                    
        # if yes, gapfill
        solution = pfba_gapfill_implementation(gf_model, gf_universal, 'biomass')
        if len(solution) >= 1:
            logger.info('reactions to add')
            logger.info(len(solution))
            for rxn_id in solution:
                add_reactions_list.append(gf_universal.reactions.get_by_id(rxn_id).copy())
                gf_mod_list2.append(gf_universal.reactions.get_by_id(rxn_id).copy())
    df = pd.DataFrame({'rxns_added':gf_mod_list2})
    df.to_csv('gapfilling_additions_no_ortho_{0}_species_biomass.csv'.format(SPECIES_ID_old))
    logger.info("wrote species biomass file")

add_reactions_list2 = list(set(hf.flatten_mixed_list([gf_mod_list1,gf_mod_list2])))
logger.info('going to add these reactions for biomass gapfill:')
logger.info(add_reactions_list2)

for rxn in add_reactions_list2:
    if rxn.id not in [r.id for r in model.reactions]:
        mets = [x.metabolites for x in [rxn]]
        all_keys = set().union(*(d.keys() for d in mets))
        for key in all_keys:
            if key.id not in [m.id for m in model.metabolites]:
                model.add_metabolites([key.copy()])
        model.add_reactions([rxn.copy()])
logger.info('added reactions for biomass gapfill')
model.repair()

os.chdir(model_path)
cobra.io.save_json_model(model, 'gf_no_ortho_'+SPECIES_ID+'.json')
cobra.io.write_sbml_model(model,'gf_no_ortho_'+SPECIES_ID+'.xml')
    
    
