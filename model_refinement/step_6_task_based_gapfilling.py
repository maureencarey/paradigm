# Input ortho_speciesname.json and iPfal19.json

import cobra
import pandas as pd
import os
import subprocess
import glob
import json
from cobra import Model, Reaction, Metabolite
from cobra.flux_analysis.parsimonious import add_pfba
import helper_functions_2 as hf2
import argparse
import logging
from datetime import datetime
from itertools import chain
from cobra.util.solver import linear_reaction_coefficients

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
logging.basicConfig(filename='step6_{}_{}.log'.format(SPECIES_ID_old,day), level=logging.INFO, filemode='w')
logger = logging.getLogger(__name__)
logger.info('BEGIN STEP 6')
    
os.chdir(model_path)
model = cobra.io.load_json_model(model_fname)
logger.info('loaded model')

# write function here since it uses logger
# pFBA based gapfilling, implementation from Greg Medlock
def pfba_gapfill_implementation(input_model, universal_model_ex, objective_reaction_id):
    # objective_reaction is a reaction id

    universal_model_pfba = universal_model_ex.copy()
    for rxn in input_model.reactions:
        if rxn.id in [r.id for r in universal_model_pfba.reactions]:
            if rxn.upper_bound > universal_model_pfba.reactions.get_by_id(rxn.id).upper_bound:
                universal_model_pfba.reactions.get_by_id(rxn.id).upper_bound = rxn.upper_bound
            if rxn.lower_bound < universal_model_pfba.reactions.get_by_id(rxn.id).lower_bound:
       	       	universal_model_pfba.reactions.get_by_id(rxn.id).lower_bound = rxn.lower_bound
        else: universal_model_pfba.add_reactions([rxn])

    # test if gapfilling is possible
    universal_model_pfba.objective = objective_reaction_id
    universal_model_pfba_test = universal_model_pfba.copy()

    # if gapfilling is not possible, don't even try. identify biomass precursors that can't be made if biomass is the objective
    sol_temp = universal_model_pfba.optimize()
    logger.info({sol_temp.status:sol_temp.objective_value})
    if sol_temp.status == 'infeasible' or sol_temp.objective_value < 0.1: # TEST IF WE EVER GO IN THIS LOOP
        logger.info('INFEASIBLE: this gapfilling problem will have no solution. demands or sinks would be needed.') 
        add_reactions_to_model = []
        if 'biomass' in objective_reaction_id:
            bm_rxn = universal_model_pfba.reactions.get_by_id(objective_reaction_id)
            total = 0
            cant = 0
            for met in bm_rxn.reactants:
                total = total+1
                if 'DM_'+met.id not in [r.id for r in universal_model_pfba.reactions]:
                    universal_model_pfba.add_boundary(met, type = "demand")
                universal_model_pfba.objective = 'DM_'+met.id
                if universal_model_pfba.slim_optimize() < 0.01: 
                    logger.info('{} cannot be synthesized for biomass'.format(met.id))
                    cant = cant+1
            logger.info('{} out of {} biomass precursors cannot be synthesized'.format(cant,total))
    else:   # if gapfilling is possible, try it  
        add_pfba(universal_model_pfba, objective = objective_reaction_id, fraction_of_optimum=0.00001)
        coef = universal_model_pfba.objective.get_linear_coefficients(universal_model_pfba.variables)
 
        # penalize adding demand and sink reactions
        dm_sk_reaction_ids = [rxn.id for rxn in universal_model_pfba.reactions if rxn.id.startswith('DM_') or rxn.id.startswith('SK_')]
        reaction_variables = (((universal_model_pfba.reactions.get_by_id(reaction).forward_variable),\
                             (universal_model_pfba.reactions.get_by_id(reaction).reverse_variable)) for reaction in dm_sk_reaction_ids)
        variables = chain(*reaction_variables)
        for variable in variables: coef[variable] = 100.0
       	
        # penalize adding exchange reactions
       	ex_reaction_ids = [rxn.id for rxn in universal_model_pfba.reactions if rxn.id.startswith('EX_')]
       	reaction_variables = (((universal_model_pfba.reactions.get_by_id(reaction).forward_variable),\
                           (universal_model_pfba.reactions.get_by_id(reaction).reverse_variable)) for reaction in ex_reaction_ids)
        variables = chain(*reaction_variables)
        for variable in variables: coef[variable] = 10.0   
        
        # remove penalty for reactions in original model # CHECK IF ZEROS ARE EXCLUDED BY DEFAULT
        original_reaction_ids = [rxn.id for rxn in input_model.reactions]
        reaction_variables = (((universal_model_pfba.reactions.get_by_id(reaction).forward_variable),\
                           (universal_model_pfba.reactions.get_by_id(reaction).reverse_variable)) for reaction in original_reaction_ids)
        variables = chain(*reaction_variables)
        for variable in variables: coef[variable] = 0.1
        universal_model_pfba.objective.set_linear_coefficients(coef)

        solution = universal_model_pfba.optimize()
        if solution.status == 'infeasible':
            logger.info('pFBA gapfilling for {} is infeasible!'.format(objective_reaction_id))
        else: 
            logger.info('pFBA solution is: {}'.format(solution.status))
        get_fluxes = set([r.id for r in universal_model_pfba.reactions]) - set([rxn.id for rxn in input_model.reactions])
        add_reactions_to_model = [rxn for rxn in get_fluxes if abs(solution.fluxes[rxn]) > 1E-10]

        test_model = input_model.copy()
        logger.info('add these (rxn ids): {}'.format(add_reactions_to_model))
        add_reactions_list = [universal_model_pfba_test.reactions.get_by_id(r).copy() for r in add_reactions_to_model]
        test_model.add_reactions(add_reactions_list)
        test_model.objective = objective_reaction_id
        sol = test_model.optimize()
        if sol.status == 'infeasible':
            logger.info('double checking pFBA solution: pFBA solution is still infeasible after double checking!')
            add_reactions_to_model = []
        else:
            logger.info('double checking pFBA solution: solution is {}'.format(sol.objective_value))
            if sol.objective_value > 0.0: 
                logger.info('double checking pFBA solution: solution is ok: {}'.format(sol.objective_value))
            else: 
                add_reactions_to_model = []
                logger.info('double checking pFBA solution: INFEASIBLE: pFBA failed to find a solution')
    return(add_reactions_to_model)

os.chdir(model_path)
universal_model = cobra.io.read_sbml_model("extended_universal_model_for_gapfilling.xml")
iPfal19	= cobra.io.read_sbml_model("iPfal19.xml")

#logger.info([r.id for r in iPfal19.reactions if r.id not in [rxn.id for rxn in universal_model.reactions]])

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

if SPECIES_ID in plasmodb: # Plasmodium = cytosol, extracellular, mitochondrdia, apicoplast, food vacuole
    compartment = ["_c","_e","_m","_ap","_fv"]
elif SPECIES_ID in tritrypdb: # Leishmania = cytosol, extracellular, mitochondrdia, kinetoplast, glycosome
    compartment = ["_c","_e","_m","_k","_glc"]
elif SPECIES_ID in cryptodb: # Cryptosporidium = cytosol, extracellular, pseudomitochondria (USE MITO)
    compartment = ["_c","_e","_m"]
elif SPECIES_ID in toxodb: # Toxoplasma = cytosol, extracellular, mitochondrdia, apicoplast
    compartment = ["_c","_e","_ap","_m"]
elif SPECIES_ID in giardiadb: # Giardia, Entamoeba = cytosol, extracellular
    compartment = ["_c","_e"]
elif SPECIES_ID in amoebadb:
    compartment = ["_c","_e"]
elif SPECIES_ID in microsporidiadb: # I haven't researched this
    compartment = ["_c","_e"]
elif SPECIES_ID in piroplasmadb: # I haven't researched this
    compartment = ["_c","_e"]
elif SPECIES_ID in trichdb: # I haven't researched this
    compartment = ["_c","_e"]
else:
    logger.info('ERROR: {} is not in any database list'.format(SPECIES_ID))
    compartment = ["_c","_e"]

logger.info('reading experimental data')
# experimental data
os.chdir(data_path)
genome_ids = pd.read_csv("auxotrophies_mapping_to_genomeID.csv",header = None).T
new_header = genome_ids.iloc[0].copy()
genome_ids = genome_ids[1:].copy()
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
       	rxn = iPfal19.reactions.get_by_id(rxn_id).copy()
        universal_model.add_reactions([rxn])
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

logger.info([r.id for r in iPfal19.reactions if r.id not in [rxn.id for rxn in universal_model.reactions]])

mets_to_prod = list()
mets_to_consume = list()
if sum(genome_ids['strain_ID'].isin([SPECIES_ID])): # does the model have any tasks?
    idx = genome_ids.index[genome_ids['strain_ID'] == SPECIES_ID].tolist()
    sp = genome_ids.loc[idx]['species'] # get species name
    species_string = sp[sp.index[0]]
    
    # tasks for that species
    gf_species = gapfilling_tasks[["Metabolite",species_string]].copy()
    indx = gf_species[species_string].notnull()
    if sum(indx) >0:
        gf_species.loc[indx]
        
        # types of tasks
        consumption_species = gf_species[gf_species[species_string].str.contains(\
        "uptake|essential|both|rescue")==True]
        production_species = gf_species[gf_species[species_string].str.contains(\
        "prod|both")==True]
        
        mets_to_prod = production_species['Metabolite'].tolist()
        mets_to_consume = consumption_species['Metabolite'].tolist()

model.solver = 'glpk'
if 'biomass' in [r.id for r in model.reactions]:
    model.objective = 'biomass' # never using this, just a test
    logger.info('was able to set objective to biomass')
elif 'generic_biomass' in [r.id for r in model.reactions]:
    model.objective = 'generic_biomass'
    logger.info('was able to set objective to generic biomass')
add_reactions_list = list()
add_mets_list = list()

#logger.info('TESTING')
#bm = iPfal19.reactions.biomass.copy()
#universal_model.add_reactions([bm])
#universal_model.objective = 'biomass'
#logger.info(universal_model.slim_optimize())
#universal_model.remove_reactions([universal_model.reactions.biomass])
#logger.info(len(universal_model.reactions))

# remove reactions that are not in eligible compartments
universal_model_for_species = universal_model.copy()
for rxn in universal_model_for_species.reactions:
    rxn_metabolite_list = [met.id for met in rxn.metabolites]
    # CHECK set(metabolite.compartment) # check if attribute is not blank
    rxn_metabolite_comp_list = [hf2.get_comp(universal_model_for_species,m) for m in rxn_metabolite_list]
    if len(hf2.unaccept_comp_intersection(rxn_metabolite_comp_list, compartment))>0:
        universal_model_for_species.remove_reactions([rxn])

#logger.info([r.id for r in iPfal19.reactions if r.id not in [rxn.id for rxn in universal_model.reactions]])

#logger.info('TESTING')
#bm = iPfal19.reactions.biomass.copy()
#universal_model_for_species.add_reactions([bm])
#universal_model_for_species.objective = 'biomass'
#logger.info(universal_model_for_species.slim_optimize())
#universal_model_for_species.remove_reactions([universal_model_for_species.reactions.biomass])

produce_met = list()
consume_met = list()
all_mets = list()
cannot_gapfill = list()
if mets_to_prod:
    for met in mets_to_prod:
        if met+'_c' in universal_model_for_species.metabolites:
            produce_met.append(met)
            all_mets.append(met)
        else:
            cannot_gapfill.append(met)
if mets_to_consume:
    for met in mets_to_consume:
        if met+'_c' in universal_model_for_species.metabolites:
            consume_met.append(met)
            all_mets.append(met)
        else:
            cannot_gapfill.append(met)
all_mets = list(set(all_mets))

logger.info('all mets = {}'.format(all_mets))

if len(all_mets) > 0 : logger.info('there are mets to gapfill for')
else: logger.info('there are no mets to gapfill for')

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
            gf_model.add_metabolites([gf_universal.metabolites.get_by_id(extracel_met).copy()])
            add_mets_list.append(gf_universal.metabolites.get_by_id(extracel_met).copy())
        else:
            met = Metabolite(extracel_met,
                            formula = gf_universal.metabolites.get_by_id(cyto_met).formula,
                            name = gf_universal.metabolites.get_by_id(cyto_met).name,
                            compartment='e')
            gf_model.add_metabolites([met])
            add_mets_list.append(met)

    # add objective reaction
    if 'DM_'+cyto_met not in [r.id for r in gf_model.reactions]:
        gf_model.add_boundary(gf_model.metabolites.get_by_id(cyto_met), type = "demand")
    gf_model.objective = 'DM_'+cyto_met

    # add exchange
    if 'EX_'+extracel_met not in gf_model.reactions:
        gf_model.add_boundary(gf_model.metabolites.get_by_id(extracel_met),type="exchange")
        add_reactions_list.append(gf_model.reactions.get_by_id('EX_'+extracel_met).copy())

    if met in list(produce_met):
        gf_model.reactions.get_by_id('EX_'+extracel_met).lower_bound = 0.
        gf_model.reactions.get_by_id('EX_'+extracel_met).upper_bound = 0.
        if not (gf_model.slim_optimize() > 0.01): # gapfill if can't produce
            solution = pfba_gapfill_implementation(gf_model, gf_universal,'DM_'+cyto_met)
            logger.info('solution to produce:')
            logger.info(solution)
            if len(solution) > 1:
                for rxn_id in solution: # should never have DM rxn here, can streamline, double check
                    add_reactions_list.append(gf_universal.reactions.get_by_id(rxn_id).copy())
        else:
            logger.info('no gapfilling needed')

    if met in list(consume_met):
        gf_model.reactions.get_by_id('EX_'+extracel_met).lb = -1000.
        gf_model.reactions.get_by_id('EX_'+extracel_met).ub = 1000.
        t = False
        if not (gf_model.slim_optimize() >= 0.01): # optimize and no growth
            # add import reaction (don't use gapfill because we want an import,
            # not biosynthetic solution)
            if extracel_met in gf_universal.metabolites:
                extracellular_reactions = gf_universal.metabolites.get_by_id(extracel_met).reactions
                intracellular_reactions = gf_universal.metabolites.get_by_id(cyto_met).reactions
                possible_rxns = list(extracellular_reactions.intersection(intracellular_reactions))
                if len(possible_rxns) >0:
                    if len(possible_rxns) > 1:
                        rxn_keep = ''
                        for rxn in possible_rxns: # add reaction that has no other substrates
                            if len(rxn.metabolites) == 2: rxn_keep = rxn
                            else: pass
                        if rxn_keep != '': add_reactions_list.append(rxn_keep)
                        else: add_reactions_list.append(possible_rxns[0])
                        logger.info('solution to consume')
                        logger.info(rxn_keep)
                    else:
                        add_reactions_list.append(possible_rxns[0])
                        logger.info('solution to consume:')
                        logger.info(possible_rxns[0])
                else:
                    t = True # if no transport in universal, make transport rxn
            else:
                t = True # if extracel met isnt in universal, make transport rxn
        else:
            logger.info('no gapfilling needed')

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

add_reactions_list = hf2.flatten_mixed_list(add_reactions_list)
df = pd.DataFrame({'rxns_added':add_reactions_list})
df.to_csv('gapfilling_additions_{0}_tasks.csv'.format(SPECIES_ID))

for met in add_mets_list:
    if met.id not in [m.id for m in model.metabolites]:
        model.add_metabolites([met])
for rxn in add_reactions_list:
    if rxn.id not in [r.id for r in model.reactions]:
        model.add_reactions([rxn])

logger.info("---------------------------------------------------")
logger.info('no. reactions added:')
logger.info(len(set(add_reactions_list)))
logger.info('no. metabolites added:')
logger.info(len(set(add_mets_list)))

#model.solver = 'glpk'
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

    logger.info([r.id for r in iPfal19.reactions if r.id not in [rxn.id for rxn in gf_universal.reactions]])

    logger.info('TESTING')
    bm = iPfal19.reactions.biomass.copy()
    gf_universal.add_reactions([bm])
    gf_universal.objective = 'biomass'
    logger.info(gf_universal.slim_optimize())
    gf_universal.remove_reactions([gf_universal.reactions.biomass])

    gf_model.objective = 'generic_biomass'
    if not (gf_model.slim_optimize() > 0.01): # gapfill if can't produce
        solution = pfba_gapfill_implementation(gf_model, gf_universal, 'generic_biomass')
        if len(solution) >= 1:
            for rxn_id in solution:
                add_reactions_list.append(gf_universal.reactions.get_by_id(rxn_id).copy())
                gf_mod_list1.append(gf_universal.reactions.get_by_id(rxn_id).copy())
    df = pd.DataFrame({'rxns_added':gf_mod_list1})
    df.to_csv('gapfilling_additions_{0}_generic_biomass.csv'.format(SPECIES_ID_old))
    logger.info("wrote generic biomass file")
else:
    logger.info('error: no generic biomass reaction')

if 'biomass' in [r.id for r in model.reactions]:
    gf_model = model.copy()
    gf_model.reactions.get_by_id('biomass').lb = 0.
    gf_model.reactions.get_by_id('biomass').ub = 1000.
    if 'generic_biomass' in [r.id for r in gf_model.reactions]:
        gf_model.remove_reactions([gf_model.reactions.get_by_id('generic_biomass')])
    gf_universal = universal_model_for_species.copy()

    logger.info([r.id for r in iPfal19.reactions if r.id not in [rxn.id for rxn in gf_universal.reactions]])

    logger.info('TESTING')
    bm = iPfal19.reactions.biomass.copy()
    gf_universal.add_reactions([bm])
    gf_universal.objective = 'biomass'
    logger.info(gf_universal.slim_optimize())
    gf_universal.remove_reactions([gf_universal.reactions.biomass])

    gf_model.objective = 'biomass'
    if not (gf_model.slim_optimize() > 0.01): # gapfill if can't produce
        solution = pfba_gapfill_implementation(gf_model, gf_universal, 'biomass')
        if len(solution) >= 1:
            for rxn_id in solution:
                add_reactions_list.append(gf_universal.reactions.get_by_id(rxn_id).copy())
                gf_mod_list2.append(gf_universal.reactions.get_by_id(rxn_id).copy())
    df = pd.DataFrame({'rxns_added':gf_mod_list2})
    df.to_csv('gapfilling_additions_{0}_species_biomass.csv'.format(SPECIES_ID_old))
    if SPECIES_ID in plasmodb:
        logger.info("wrote species biomass file")
    else: logger.info("ERROR: wrote the plasmodium species biomass file but this is not a plasmodium model")
else:
    if SPECIES_ID in plasmodb:
        logger.info('error: no plasmodium biomass reaction')
    else: logger.info('no plasmodium biomass reaction, which is what was intended')

add_reactions_list2 = list(set(hf2.flatten_mixed_list([gf_mod_list1,gf_mod_list2])))
logger.info('going to add these reactions for biomass gapfill:')
logger.info(add_reactions_list2)
for rxn in add_reactions_list2:
    if rxn.id not in [r.id for r in model.reactions]:
        model.add_reactions([rxn])
logger.info('added reactions for biomass gapfill')

os.chdir(model_path)
cobra.io.save_json_model(model, 'gf_'+SPECIES_ID+'.json')
print(hippo)
cobra.io.write_sbml_model(model,'gf_'+SPECIES_ID+'.xml')
    
    
