# Input with_biomass_speciesname.json and iPfal19.json

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

data_path = "/home/mac9jc/paradigm/data"
model_path = "/home/mac9jc/paradigm/models"
os.chdir(model_path)

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

def intersection(rxn_list_compartments, acceptable_compartments):
    temp = set(acceptable_compartments)
    unacceptable_comp = [value for value in rxn_list_compartments if value not in temp]
    return(unacceptable_comp)

def get_comp(model,met_id):

    # get compartment associated with a metabolite(s)
    if met_id in [met.id for met in model.metabolites]:
        for m in [model.metabolites.get_by_id(met_id)]:
            if m.id.endswith('_c') or m.id.endswith('_e') or m.id.endswith('_f') or \
            m.id.endswith('_g') or m.id.endswith('_h') or m.id.endswith('_i') or \
            m.id.endswith('_l') or m.id.endswith('_m') or m.id.endswith('_n') or \
            m.id.endswith('_p') or m.id.endswith('_r') or m.id.endswith('_s') or \
            m.id.endswith('_u') or m.id.endswith('_v') or m.id.endswith('_x'):
                id_withou_c = m.id[-2:]
            elif m.id.endswith('_cx') or m.id.endswith('_um') or m.id.endswith('_im') or \
             m.id.endswith('_cm') or m.id.endswith('_ap') or m.id.endswith('_fv'):
                id_withou_c = m.id[-3:]
            else:
                print('unknown compartment')
                print(m.id)
                id_withou_c = ''
    else:
        id_withou_c = ''
    return(id_withou_c)
    
# modified for Rivanna: read in the models
model_dict = {}
model_dict[SPECIES_ID] = cobra.io.load_json_model(model_fname)
logger.info('loaded model')

for species, model in model_dict.items():
    if 'biomass' not in [r.id for r in model.reactions]:
        logger.info('biomass not in reactions anymore')

os.chdir(model_path)
universal_model = cobra.io.load_json_model('universal_model_updated.json')

# extend universal by curated model
pf_model = cobra.io.load_json_model('iPfal19_updated.json')
len_univ_rxns = len(universal_model.reactions)
for rxn in pf_model.reactions:
    if rxn.id not in [r.id for r in universal_model.reactions]:
        if len(set(['hb_c','hb_e','hemozoin_c','hemozoin_e','hemozoin_fv']).intersection(set([met.id for met in rxn.metabolites.keys()]))) == 0:
            universal_model.add_reactions([rxn.copy()]) # extend universal by Pf reactions, but remove gene IDs
            if len(rxn.genes) > 0:
                genes = rxn.genes
                universal_model.reactions.get_by_id(rxn.id).gene_reaction_rule = ''
                for gene in genes:
                    if len(universal_model.genes.get_by_id(gene.id).reactions) == 0:
                           gene_list = [universal_model.genes.get_by_id(gene.id)]
                           cobra.manipulation.delete.remove_genes(universal_model, gene_list, remove_reactions=True)
    else: #rxn.id IS in [r.id for r in universal_model but in PF its reversible and in universal its irreversible, make reversible
        if rxn.lower_bound < universal_model.reactions.get_by_id(rxn.id).lower_bound:
            universal_model.reactions.get_by_id(rxn.id).lower_bound = rxn.lower_bound
        if rxn.upper_bound > universal_model.reactions.get_by_id(rxn.id).upper_bound:
            universal_model.reactions.get_by_id(rxn.id).upper_bound = rxn.upper_bound
if len(universal_model.reactions) <= len_univ_rxns:
    logger.info('ERROR - universal model does not have Pf reactions added!')

## prep for removing universal reactions that are in the wrong compartment for this model
# database mapping - this is for release 42, must update if using a different EuPathDB release
plasmodb = ["PadleriG01","PbergheiANKA","PbillcollinsiG01","PblacklockiG01","Pchabaudichabaudi","PcoatneyiHackeri","PcynomolgiB","PcynomolgiM","Pfalciparum3D7","Pfalciparum7G8","PfalciparumCD01","PfalciparumDd2","PfalciparumGA01","PfalciparumGB4","PfalciparumGN01","PfalciparumHB3","PfalciparumIT","PfalciparumKE01","PfalciparumKH01","PfalciparumKH02","PfalciparumML01","PfalciparumSD01","PfalciparumSN01","PfalciparumTG01","PfragileNilgiri","PgaboniG01","PgaboniSY75","Pgallinaceum8A","PinuiSanAntonio1","PknowlesiH","PknowlesiMalayanPk1A","PmalariaeUG01","PovalecurtisiGH01","PpraefalciparumG01","PreichenowiCDC","PreichenowiG01","PrelictumSGS1-like","PvinckeipetteriCR","Pvinckeivinckeivinckei","Pvivax-likePvl01","PvivaxP01","PvivaxSal1","Pyoeliiyoelii17X","Pyoeliiyoelii17XNL","PyoeliiyoeliiYM"]
cryptodb = ["Candersoni30847","Chominis30976","ChominisTU502","ChominisTU502_2012","ChominisUdeA01","CmeleagridisUKMEL1","CmurisRN66","CparvumIowaII","CtyzzeriUGA55", "Cubiquitum39726","CveliaCCMP2878", "GniphandrodesUnknown", "VbrassicaformisCCMP3155"]
giardiadb = ["GintestinalisAssemblageADH", "GintestinalisAssemblageAWB", "GintestinalisAssemblageBGS", "GintestinalisAssemblageBGS_B", "GintestinalisAssemblageEP15", "SsalmonicidaATCC50377"]
tritrypdb = ["BayalaiB08-376","BsaltansLakeKonstanz","CfasciculataCfCl","EmonterogeiiLV88","LaethiopicaL147", "LamazonensisMHOMBR71973M2269","LarabicaLEM1108", "LbraziliensisMHOMBR75M2903", "LbraziliensisMHOMBR75M2904", "LdonovaniBPK282A1","LdonovaniCL-SL", "LenriettiiLEM3045", "LgerbilliLEM452","LinfantumJPCM5", "LmajorFriedlin", "LmajorLV39c5", "LmajorSD75.1","LmajorSD75", "LmexicanaMHOMGT2001U1103", "LpanamensisMHOMCOL81L13","LpanamensisMHOMPA94PSC1", "LpyrrhocorisH10", "LseymouriATCC30220", "LspMARLEM2494", "LtarentolaeParrotTarII", "LtropicaL590", "LturanicaLEM423", "PconfusumCUL13","TbruceigambienseDAL972", "TbruceiLister427", "TbruceiLister427_2018", "TbruceiTREU927", "TcongolenseIL3000", "TcruziCLBrener", "TcruziCLBrenerEsmeraldo-like", "TcruziCLBrenerNon-Esmeraldo-like", "TcruziDm28c2014","TcruziDm28c2017","TcruziDm28c2018", "TcruzimarinkelleiB7", "TcruziSylvioX10-1", "TcruziSylvioX10-1-2012","TcruziTCC","TevansiSTIB805", "TgrayiANR4", "TrangeliSC58", "TvivaxY486", "TtheileriEdinburgh"]
trichdb = ["TvaginalisG3"]
amoebadb = ["AcastellaniiNeff", "EdisparSAW760", "EhistolyticaHM1IMSS-A", "EhistolyticaHM1IMSS-B", "EhistolyticaHM1IMSS", "EhistolyticaHM3IMSS", "EhistolyticaKU27", "EinvadensIP1", "EmoshkovskiiLaredo", "EnuttalliP19", "NfowleriATCC30863"]
toxodb = ["CcayetanensisCHN_HEN01", "CsuisWienI","EacervulinaHoughton", "EbrunettiHoughton", "EfalciformisBayerHaberkorn1970", "EmaximaWeybridge", "EmitisHoughton", "EnecatrixHoughton", "EpraecoxHoughton", "EtenellaHoughton", "HhammondiHH34", "NcaninumLIV", "SneuronaSN3", "SneuronaSOSN1", "TgondiiARI", "TgondiiFOU", "TgondiiGAB2-2007-GAL-DOM2", "TgondiiGT1", "TgondiiMAS", "TgondiiME49", "Tgondiip89", "TgondiiRH", "TgondiiRUB", "TgondiiTgCatPRC2", "TgondiiVAND", "TgondiiVEG"]
microsporidiadb = ["AalgeraePRA109", "AalgeraePRA339", "AspWSBS2006","EaedisUSNM41457", "EbieneusiH348", "EcanceriGB1","EcuniculiEC1", "EcuniculiEC2", "EcuniculiEC3","EcuniculiEcunIII-L","EcuniculiGBM1", "EhellemATCC50504", "EhellemSwiss", "EhepatopenaeiTH1","EintestinalisATCC50506", "EromaleaeSJ2008","Heriocheircanceri","HeriocheirGB1", "MdaphniaeUGP3", "NausubeliERTm2", "NausubeliERTm6", "NbombycisCQ1", "NceranaeBRL01","NceranaePA08_1199","NdisplodereJUm2807","NparisiiERTm1", "NparisiiERTm3", "OcolligataOC4", "PneurophiliaMK1", "Slophii42_110", "ThominisUnknown", "VcorneaeATCC50505", "Vculicisfloridensis"]
piroplasmadb = ["BbigeminaBOND", "BbovisT2Bo", "Bdivergens1802A","BmicrotiRI","BovataMiyake", "CfelisWinnie", "TannulataAnkara", "TequiWA", "TorientalisShintoku", "TparvaMuguga"]

compartment_dictionary = dict()
for species in model_dict.keys():
    # Plasmodium = cytosol, extracellular, mitochondrdia, apicoplast, food vacuole
    if species in plasmodb:
        model_compartments = ["_c","_e","_m","_ap","_fv"]
    # Leishmania = cytosol, extracellular, mitochondrdia, kinetoplast, glycosome
    elif species in tritrypdb:
        model_compartments = ["_c","_e","_m","_k","_glc"]
    # Cryptosporidium = cytosol, extracellular, pseudomitochondria (USE MITO)
    elif species in cryptodb:
        model_compartments = ["_c","_e","_m"]
    # Toxoplasma = cytosol, extracellular, mitochondrdia, apicoplast
    elif species in toxodb:
        model_compartments = ["_c","_e","_ap","_m"]
    # Giardia, Entamoeba = cytosol, extracellular
    elif species in giardiadb or species in amoebadb:
        model_compartments = ["_c","_e"]
    else:
        model_compartments = ["_c","_e"]
    compartment_dictionary[species] = model_compartments

logger.info('loaded universal')
os.chdir(data_path)
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
    print(solution)
    print(dir(solution))
    add_reactions_to_model = [rxn for rxn in get_fluxes if abs(solution.fluxes[rxn]) > 1E-8]

    # double check
    logging.info('double checking pFBA solution')
    test_model = input_model.copy()
    add_reactions_list = [universal_model_pfba.reactions.get_by_id(r).copy() for r in add_reactions_to_model]
    test_model.add_reactions(add_reactions_list)
    sol = test_model.optimize()
    if sol.status == 'infeasible':
        logging.info('pFBA solution is infeasible!')
    return(add_reactions_to_model)

GENERIC_BIOMASS = 'generic_biomass'
P_BIOMASS = 'biomass'

for species, model in model_dict.items():
    
    logger.info(species)
    model.solver = 'glpk'
    if 'biomass' in [r.id for r in model.reactions]:
        model.objective = 'biomass'
        logger.info('was able to set objective')
    add_reactions_list = list()
    add_mets_list = list()

    # remove reactions that are not in eligible compartments
    compartments = compartment_dictionary[species]
    universal_model_for_species = universal_model.copy()
    for rxn in universal_model_for_species.reactions:
        # rxn_metabolite_list = '\t'.join([m.id for m in rxn.metabolites])
        rxn_metabolite_list = [get_comp(universal_model_for_species,m.id) for m in rxn.metabolites]
        if len(intersection(rxn_metabolite_list, compartments))>0:
            universal_model_for_species.remove_reactions([rxn])    

    if species in tasks_dict.keys():

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

        for met in all_mets:
            logger.info('')
            logger.info(met)
            extracel_met = met+'_e'
            cyto_met = met+'_c'
    
            gf_model = model.copy()
            gf_universal = universal_model_for_species.copy()
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

    # is it used in any reaction intracellularly? 
    # if not, cannot do anything about it, except curate in future

        add_reactions_list = flatten_mixed_list(add_reactions_list)
        df = pd.DataFrame({'rxns_added':add_reactions_list})
        df.to_csv('gapfilling_additions_no_ortho_{0}_tasks.csv'.format(species))

        for met in add_mets_list:
            if met.id not in [m.id for m in model.metabolites]:
                model.add_metabolites([met])
        for rxn in add_reactions_list:
            if rxn.id not in [r.id for r in model.reactions]:
                model.add_reactions([rxn])
    else:
        logger.info("no data to gapfill for")

    logger.info("---------------------------------------------------")
    logger.info(species)
    logger.info('no. reactions added:')
    logger.info(len(set(add_reactions_list)))
    logger.info('no. metabolites added:')
    logger.info(len(set(add_mets_list)))
    # ADD FIRST
    model_dict[species] = model
    # save models

    logger.info(species)
    model.solver = 'glpk'
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
        gf_universal = universal_model_for_species.copy()
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
        df.to_csv('gapfilling_additions_no_ortho_{0}_generic_biomass.csv'.format(SPECIES_ID_old))
        logger.info("wrote generic bioamss file")
    else:
        logger.info('error: no generic biomass reaction')
    
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
        df.to_csv('gapfilling_additions_no_ortho_{0}_species_biomass.csv'.format(SPECIES_ID_old))
        logger.info("wrote species biomass file")

    add_reactions_list2 = list(set(flatten_mixed_list([gf_mod_list1,gf_mod_list2])))
    logger.info('going to add these reactions:')
    logger.info(add_reactions_list2)
    for rxn in add_reactions_list2:
        if rxn.id not in [r.id for r in model.reactions]:
            model.add_reactions([rxn])
    logger.info('added reactions')

    os.chdir(model_path)
    cobra.io.save_json_model(model, 'gf_no_ortho_'+SPECIES_ID+'.json')
    cobra.io.write_sbml_model(model,'gf_no_ortho_'+SPECIES_ID+'.xml')
    
    
