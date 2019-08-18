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

# write function here since it useses logger
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
    if universal_model_pfba.slim_optimize() < 0.1: 
        logger.info('INFEASIBLE: this gapfilling problem will have no solution') 
        skip = 1
    else: skip = 0

    add_pfba(universal_model_pfba)

    # penalize adding demand reactions
    coef = universal_model_pfba.objective.get_linear_coefficients(universal_model_pfba.variables)
    for key,value in coef.items():
        if key.name.startswith('DM_') or key.name.startswith('SK_'):
            coef[key] = 100.
        elif key.name.startswith('EX_'):
            coef[key] = 10.
    for x in universal_model_pfba.variables:
        if x.name.startswith('DM_') or x.name.startswith('SK_') or x.name.startswith('EX_'):
            universal_model_pfba.objective.set_linear_coefficients(coef)

    # remove penalty for reactionsn in original model
    coef2 = universal_model_pfba.objective.get_linear_coefficients(universal_model_pfba.variables)
    for key,value in coef2.items():
        if key.name in [rxn.id for rxn in input_model.reactions]:
            coef2[key] = 0.
    for x in universal_model_pfba.variables:
        if x.name in [rxn.id for rxn in input_model.reactions]:
            universal_model_pfba.objective.set_linear_coefficients(coef2)

    universal_model_pfba.reactions.get_by_id(objective_reaction_id).lower_bound = 0.01

    solution = universal_model_pfba.optimize()
    if solution.status == 'infeasible' and skip == 0:
        logger.info('pFBA gapfilling for {} is infeasible!'.format(objective_reaction_id))
    get_fluxes = set([r.id for r in universal_model_pfba.reactions]) - set([rxn.id for rxn in input_model.reactions])
    add_reactions_to_model = [rxn for rxn in get_fluxes if abs(solution.fluxes[rxn]) > 1E-8]

    # double check
    logger.info('double checking pFBA solution')
    test_model = input_model.copy()
    add_reactions_list = [universal_model_pfba.reactions.get_by_id(r).copy() for r in add_reactions_to_model]
    test_model.add_reactions(add_reactions_list)
    sol = test_model.optimize()
    if sol.status == 'infeasible':
        logger.info('pFBA solution is still infeasible after double checking!')
    else:
       	f = test_model.slim_optimize()
       	if f > 0.0 and skip == 0:
            logger.info('solution is ok: {}'.format(test_model.slim_optimize()))
       	#else:
       	elif skip == 0:
            logger.info('INFEASIBLE: pFBA failed to find a solution')
    if skip == 1: add_reactions_to_model = []
    return(add_reactions_to_model)

if 'biomass' not in [r.id for r in model.reactions]:
    if 'generic_biomass' not in [r.id for r in model.reactions]:
        logger.info('biomass not in reactions anymore')

os.chdir(model_path)
universal_model = cobra.io.read_sbml_model("extended_universal_model_for_gapfilling.xml")
iPfal19	= cobra.io.read_sbml_model("iPfal19.xml")

logger.info([r.id for r in iPfal19.reactions if r.id not in [rxn.id for rxn in universal_model.reactions]])

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
else:
    compartment = ["_c","_e"]

logger.info(compartment)
logger.info(len(hf2.intersection(['sample_c','sample_e','sample_p'], compartment))>0)
logger.info(hf2.intersection(['sample_c','sample_e','sample_p'], compartment))
logger.info([hf2.get_comp(universal_model,m.id) for m in universal_model.reactions[0].metabolites])
logger.info(t)

# experimental data
os.chdir(data_path)
genome_ids = pd.read_csv("auxotrophies_mapping_to_genomeID.csv",header = None).T
new_header = genome_ids.iloc[0]
genome_ids = genome_ids[1:]
genome_ids.columns = new_header
met_ids = pd.read_csv("auxotrophies_mapping_to_metID.csv")
gapfilling_tasks = pd.read_excel("auxotrophies_references.xlsx",skiprows = 1)
idx2 = met_ids['BiGG'].notnull()
# temp_met_ids = met_ids.loc[idx2]
temp_met_ids = met_ids.loc[met_ids['BiGG'].notnull()].copy()
# met_dict = pd.Series(temp_met_ids.BiGG.values, index = temp_met_ids.Metabolite).to_dict()
met_dict = pd.Series(met_ids.loc[idx2].BiGG.values, index = met_ids.loc[idx2].Metabolite).to_dict()

# change met ids to model met ids if possible
for x in gapfilling_tasks.index:
    if gapfilling_tasks['Metabolite'].iloc[x] in met_dict.keys():
        gapfilling_tasks['Metabolite'].iloc[x] = met_dict[gapfilling_tasks['Metabolite'].iloc[x]]
        gapfilling_tasks.to_csv('temp_gf.csv')

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

logger.info('TESTING')
bm = iPfal19.reactions.biomass.copy()
universal_model.add_reactions([bm])
universal_model.objective = 'biomass'
logger.info(universal_model.slim_optimize())
universal_model.remove_reactions([universal_model.reactions.biomass])

logger.info(compartment)
# remove reactions that are not in eligible compartments
universal_model_for_species = universal_model.copy()
for rxn in universal_model_for_species.reactions:
    rxn_metabolite_list = [hf2.get_comp(universal_model_for_species,m.id) for m in rxn.metabolites]
    if len(hf2.intersection(rxn_metabolite_list, compartment))>0:
        universal_model_for_species.remove_reactions([rxn])

for rxn in iPfal19.reactions:
    if rxn.id not in [r.id for r in universal_model_for_species.reactions]:
        print(rxn.id)


logger.info('TESTING')
bm = iPfal19.reactions.biomass.copy()
universal_model_for_species.add_reactions([bm])
universal_model_for_species.objective = 'biomass'
logger.info(universal_model_for_species.slim_optimize())
universal_model_for_species.remove_reactions([universal_model_for_species.reactions.biomass])

print(t)

#logger.info([r.id for r in iPfal19.reactions if r.id not in [rxn.id for rxn in universal_model_for_species.reactions]])

#for rxn in iPfal19.reactions:
#    if rxn.id in [r.id for r in universal_model_for_species.reactions]:
#        if rxn.reaction != universal_model_for_species.reactions.get_by_id(rxn.id).reaction:
#            logger.info(rxn.id)
#            logger.info(rxn.reaction)
#            logger.info(universal_model_for_species.reactions.get_by_id(rxn.id).reaction)
#        if rxn.upper_bound !=universal_model_for_species.reactions.get_by_id(rxn.id).upper_bound or rxn.lower_bound != universal_model_for_species.reactions.get_by_id(rxn.id).lower_bound:
#            logger.info(rxn.id)

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

    logger.info('TESTING')
    bm = iPfal19.reactions.biomass.copy()
    gf_universal.add_reactions([bm])
    gf_universal.objective = 'biomass'
    logger.info(gf_universal.slim_optimize())
    gf_universal.remove_reactions([gf_universal.reactions.biomass])

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
        logger.info('add exchange')
        gf_model.add_boundary(gf_model.metabolites.get_by_id(extracel_met),type="exchange")
        add_reactions_list.append(gf_model.reactions.get_by_id('EX_'+extracel_met).copy())

    if met in list(produce_met):
        logger.info('produce')
        gf_model.reactions.get_by_id('EX_'+extracel_met).lower_bound = 0.
        gf_model.reactions.get_by_id('EX_'+extracel_met).upper_bound = 0.
        if not (gf_model.slim_optimize() > 0.01): # gapfill if can't produce
            solution = pfba_gapfill_implementation(gf_model, gf_universal,'DM_'+cyto_met)
            logger.info('solution')
            logger.info(solution)
            if len(solution) > 1:
                for rxn_id in solution:
                    if rxn_id != 'DM_'+cyto_met:
                       add_reactions_list.append(gf_universal.reactions.get_by_id(rxn_id).copy())
                    else:
                        if solution[0] != 'DM_'+cyto_met:
                           add_reactions_list.append(gf_universal.reactions.get_by_id(solution[0]).copy())
        else:
            logger.info('no gapfilling needed')

    if met in list(consume_met):
        logger.info('consume')
        gf_model.reactions.get_by_id('EX_'+extracel_met).lb = -1000.
        gf_model.reactions.get_by_id('EX_'+extracel_met).ub = 1000.
        t = False
        if not (gf_model.slim_optimize() >= 0.01): # optimize and no growth
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
else:
    logger.info("no data to gapfill for")

logger.info("---------------------------------------------------")
logger.info('no. reactions added:')
logger.info(len(set(add_reactions_list)))
logger.info('no. metabolites added:')
logger.info(len(set(add_mets_list)))

model.solver = 'glpk'
gf_mod_list1 = list()
gf_mod_list2 = list()

logger.info('generic biomass in reactions')
logger.info('generic_biomass' in [r.id for r in model.reactions])

if 'generic_biomass' in [r.id for r in model.reactions]:
    gf_model = model.copy()
    gf_model.reactions.get_by_id('generic_biomass').lb = 0.
    gf_model.reactions.get_by_id('generic_biomass').ub = 1000.
    if 'biomass' in [r.id for r in gf_model.reactions]:
        gf_model.remove_reactions([gf_model.reactions.get_by_id('biomass')])
    gf_universal = universal_model_for_species.copy()
    gf_model.objective = 'generic_biomass'
    if not (gf_model.slim_optimize() > 0.01): # gapfill if can't produce
         solution = pfba_gapfill_implementation(gf_model, gf_universal, 'generic_biomass')
         if len(solution) >= 1:
             for rxn_id in solution:
                 add_reactions_list.append(gf_universal.reactions.get_by_id(rxn_id).copy())
                 gf_mod_list1.append(gf_universal.reactions.get_by_id(rxn_id).copy())
    os.chdir(data_path)
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

#    logger.info('TESTING')
#    bm = iPfal19.reactions.biomass.copy()
#    gf_universal.add_reactions([bm])
#    gf_universal.objective = 'biomass'
#    logger.info(gf_universal.slim_optimize())
#    gf_universal.remove_reactions([gf_universal.reactions.biomass])

    gf_model.objective = 'biomass'
    if not (gf_model.slim_optimize() > 0.01): # gapfill if can't produce
         solution = pfba_gapfill_implementation(gf_model, gf_universal, 'biomass')
         if len(solution) >= 1:
             logger.info('reactions to add')
             logger.info(len(solution))
             for rxn_id in solution:
                 add_reactions_list.append(gf_universal.reactions.get_by_id(rxn_id).copy())
                 gf_mod_list2.append(gf_universal.reactions.get_by_id(rxn_id).copy())
    os.chdir(data_path)
    df = pd.DataFrame({'rxns_added':gf_mod_list2})
    df.to_csv('gapfilling_additions_{0}_species_biomass.csv'.format(SPECIES_ID_old))
    logger.info("wrote species biomass file")

add_reactions_list2 = list(set(hf2.flatten_mixed_list([gf_mod_list1,gf_mod_list2])))
logger.info('going to add these reactions for biomass gapfill:')
logger.info(add_reactions_list2)
for rxn in add_reactions_list2:
    if rxn.id not in [r.id for r in model.reactions]:
        model.add_reactions([rxn])
logger.info('added reactions for biomass gapfill')

os.chdir(model_path)
cobra.io.save_json_model(model, 'gf_'+SPECIES_ID+'.json')
cobra.io.write_sbml_model(model,'gf_'+SPECIES_ID+'.xml')
    
    
