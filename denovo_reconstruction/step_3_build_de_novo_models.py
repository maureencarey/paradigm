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
logging.basicConfig(filename='step3_{}_{}.log'.format(SPECIES_ID,day), level=logging.INFO, filemode='w')
logger = logging.getLogger(__name__)

# begin
logging.info('BEGIN STEP 3')
logging.info(SPECIES_ID)
data_path = "/home/mac9jc/paradigm/data"
model_path = "/home/mac9jc/paradigm/models"
og_path = "/home/mac9jc/paradigm/"

# universal reaction bag for model generation
os.chdir(model_path)
universal_model = cobra.io.read_sbml_model('/home/mac9jc/paradigm/models/universal/universal_model_updated.xml')
universal_model = hf.update_universal_model(universal_model)
os.chdir(data_path)
logging.info('loaded universal')

# get KEGG ids for EuPathDB annotation addition
universal_KEGG_dict = hf.transform_universal_to_KEGG(universal_model)

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
compartment_shorthand = {'cytoplasm':'_c', 'extracellular':'_e', 'periplasm':'_p', 
                         'mitochondrion':'_m', 'peroxisome/glyoxysome':'_x', 'nucleus':'_n', 
                         'endoplasmic reticulum':'_r', 'vacuole':'_v', 'golgi apparatus':'_g', 
                         'thylakoid':'_u', 'lysosome':'_l', 'eyespot':'_s', 
                         'chloroplast':'_h', 'flagellum':'_f', 'intermembrane space of mitochondria':'_im', 
                         'thylakoid membrane':'_um', 'carboxyzome':'_cx', 'not_provided_by_bigg':'_i', 
                         'cytosolic membrane':'_cm'}
universal_dict_test = [{rxn.id:
         {'reactants':[hf.met_ids_without_comp(universal_model,m.id) for m in rxn.reactants],
          'products':[hf.met_ids_without_comp(universal_model,m.id) for m in rxn.products],
          'compartment':[compartment_shorthand[x] for x in rxn.compartments]}} 
        for rxn in universal_model.reactions]
universal_dict = dict()
for mini_dict in universal_dict_test:
    for key, value in mini_dict.items():
        universal_dict[key] = value
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
cryptodb = ["Candersoni30847","Chominis30976","ChominisTU502","ChominisTU502_2012","ChominisUdeA01","CmeleagridisUKMEL1","CmurisRN66","CparvumIowaII","CtyzzeriUGA55", "Cubiquitum39726","CveliaCCMP2878", "GniphandrodesUnknown", "VbrassicaformisCCMP3155"]
giardiadb = ["GintestinalisAssemblageADH", "GintestinalisAssemblageAWB", "GintestinalisAssemblageBGS", "GintestinalisAssemblageBGS_B", "GintestinalisAssemblageEP15", "SsalmonicidaATCC50377"]
tritrypdb = ["BayalaiB08-376","BsaltansLakeKonstanz","CfasciculataCfCl","EmonterogeiiLV88","LaethiopicaL147", "LamazonensisMHOMBR71973M2269","LarabicaLEM1108", "LbraziliensisMHOMBR75M2903", "LbraziliensisMHOMBR75M2904", "LdonovaniBPK282A1","LdonovaniCL-SL", "LenriettiiLEM3045", "LgerbilliLEM452","LinfantumJPCM5", "LmajorFriedlin", "LmajorLV39c5", "LmajorSD75.1","LmajorSD75", "LmexicanaMHOMGT2001U1103", "LpanamensisMHOMCOL81L13","LpanamensisMHOMPA94PSC1", "LpyrrhocorisH10", "LseymouriATCC30220", "LspMARLEM2494", "LtarentolaeParrotTarII", "LtropicaL590", "LturanicaLEM423", "PconfusumCUL13","TbruceigambienseDAL972", "TbruceiLister427", "TbruceiLister427_2018", "TbruceiTREU927", "TcongolenseIL3000", "TcruziCLBrener", "TcruziCLBrenerEsmeraldo-like", "TcruziCLBrenerNon-Esmeraldo-like", "TcruziDm28c2014","TcruziDm28c2017","TcruziDm28c2018", "TcruzimarinkelleiB7", "TcruziSylvioX10-1", "TcruziSylvioX10-1-2012","TcruziTCC","TevansiSTIB805", "TgrayiANR4", "TrangeliSC58", "TvivaxY486", "TtheileriEdinburgh"]
trichdb = ["TvaginalisG3"]
amoebadb = ["AcastellaniiNeff", "EdisparSAW760", "EhistolyticaHM1IMSS-A", "EhistolyticaHM1IMSS-B", "EhistolyticaHM1IMSS", "EhistolyticaHM3IMSS", "EhistolyticaKU27", "EinvadensIP1", "EmoshkovskiiLaredo", "EnuttalliP19", "NfowleriATCC30863"]
toxodb = ["CcayetanensisCHN_HEN01", "CsuisWienI","EacervulinaHoughton", "EbrunettiHoughton", "EfalciformisBayerHaberkorn1970", "EmaximaWeybridge", "EmitisHoughton", "EnecatrixHoughton", "EpraecoxHoughton", "EtenellaHoughton", "HhammondiHH34", "NcaninumLIV", "SneuronaSN3", "SneuronaSOSN1", "TgondiiARI", "TgondiiFOU", "TgondiiGAB2-2007-GAL-DOM2", "TgondiiGT1", "TgondiiMAS", "TgondiiME49", "Tgondiip89", "TgondiiRH", "TgondiiRUB", "TgondiiTgCatPRC2", "TgondiiVAND", "TgondiiVEG"]
microsporidiadb = ["AalgeraePRA109", "AalgeraePRA339", "AspWSBS2006","EaedisUSNM41457", "EbieneusiH348", "EcanceriGB1","EcuniculiEC1", "EcuniculiEC2", "EcuniculiEC3","EcuniculiEcunIII-L","EcuniculiGBM1", "EhellemATCC50504", "EhellemSwiss", "EhepatopenaeiTH1","EintestinalisATCC50506", "EromaleaeSJ2008","Heriocheircanceri","HeriocheirGB1", "MdaphniaeUGP3", "NausubeliERTm2", "NausubeliERTm6", "NbombycisCQ1", "NceranaeBRL01","NceranaePA08_1199","NdisplodereJUm2807","NparisiiERTm1", "NparisiiERTm3", "OcolligataOC4", "PneurophiliaMK1", "Slophii42_110", "ThominisUnknown", "VcorneaeATCC50505", "Vculicisfloridensis"]
piroplasmadb = ["BbigeminaBOND", "BbovisT2Bo", "Bdivergens1802A","BmicrotiRI","BovataMiyake", "CfelisWinnie", "TannulataAnkara", "TequiWA", "TorientalisShintoku", "TparvaMuguga"]


# remove duplciate reactions in mulitple compartments
if SPECIES_ID in plasmodb: # Plasmodium = cytosol, extracellular, mitochondrdia, apicoplast, food vacuole
    compartment = ["_c","_e","_m","_ap","_fv"]
    database = "PlasmoDB"
    df_use_name = 'PlasmoDB_Feb92021_MetPath.txt'
    gene_string = '-t36_1-p1'
elif SPECIES_ID in tritrypdb: # Leishmania = cytosol, extracellular, mitochondrdia, kinetoplast, glycosome
    compartment = ["_c","_e","_m","_k","_glc"]
    database = "TriTrypDB"
    df_use_name = 'TriTrypDB_Feb92021_MetPath.txt'
    gene_string = ''
elif SPECIES_ID in cryptodb: # Cryptosporidium = cytosol, extracellular, pseudomitochondria (USE MITO)
    compartment = ["_c","_e","_m"]
    database = "CryptoDB"
    df_use_name = 'CryptoDB_Feb92021_MetPath.txt'
    gene_string = ''
elif SPECIES_ID in toxodb: # Toxoplasma = cytosol, extracellular, mitochondrdia, apicoplast
    compartment = ["_c","_e","_ap","_m"]
    database = "ToxoDB"
    df_use_name = 'ToxoDB_Feb92021_MetPath.txt'
    gene_string = ''
elif SPECIES_ID in giardiadb: # Giardia, Entamoeba = cytosol, extracellular
    compartment = ["_c","_e"]
    database = "GiardiaDB"
    df_use_name = 'GiardiaDB_Feb92021_MetPath.txt'
    gene_string = ''
elif SPECIES_ID in amoebadb:
    compartment = ["_c","_e"]
    database = "AmoebaDB"
    df_use_name = 'AmoebaDB_Feb92021_MetPath.txt'
    gene_string = ''
elif SPECIES_ID in microsporidiadb: # I haven't researched this
    compartment = ["_c","_e"]
    database = "MicrosporidiaDB"
    df_use_name = 'MicrosporidiaDB_Feb92021_MetPath.txt'
    gene_string = ''
elif SPECIES_ID in piroplasmadb: # I haven't researched this
    compartment = ["_c","_e"]
    database = "PiroplasmaDB"
    df_use_name = 'PiroplasmaDB_Feb92021_MetPath.txt'
    gene_string = ''
elif SPECIES_ID in trichdb: # I haven't researched this
    compartment = ["_c","_e"]
    database = "TrichDB"
    df_use_name = 'TrichDB_Feb92021_MetPath.txt'
    gene_string = ''
else:
    logging.info("error - species ID not in database lists - will cause problem")
    compartment = ["_c","_e"]

columns = ['species','reactions_removed1','reactions_added','genes_only_on_EuPathDB']
modifications = pd.DataFrame(index = [0], columns=columns)

## get EuPathDB annotations
os.chdir("/home/mac9jc/paradigm/data/VEuPathDB_KEGG")
KEGG_DB = pd.read_csv(df_use_name, sep = "\t")
KEGG_DB = KEGG_DB.loc[KEGG_DB['Exact EC Number Match'] == 'Yes']
KEGG_DB = KEGG_DB.loc[KEGG_DB['Source'] == 'KEGG']
KEGG_DB = KEGG_DB[['Gene ID','Reaction']]


# load Diamond annotation data file
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
 
g = new_model.genes
   
# now get EuPathDB annotations
#filepath = "/home/mac9jc/paradigm/data/genomes/protein/zipped_protein/" + database +"-50_" + SPECIES_ID + "_AnnotatedProteins.fasta" 
filepath = "/home/mac9jc/paradigm/data/genomes/protein/zipped_protein/" + SPECIES_ID + "_annotatedProteins.fasta"
file1 = open(filepath, 'r')
count = 0

add_EuPath = list()
EuPath_reaction_gene_dict = dict()
while True:
    # Get next line from file
    line = file1.readline()
    # if line is empty, end of file is reached
    if not line:
        break
    # extract only gene ids, not protein sequences
    if '>' in line:
        protein = line.split(' | ')[0][1:]
        #gene_id = hf.prune_protein_to_gene_id(protein,'-t36_1-p1')
        gene_id = protein
        EuPath_reaction_gene_dict[hf.transcript_to_gene_id(gene_id)] = hf.get_KEGG_id(hf.transcript_to_gene_id(gene_id),KEGG_DB)
        for KID in EuPath_reaction_gene_dict[hf.transcript_to_gene_id(gene_id)]:
            RXNs = hf.get_rxn_from_KEGG(KID, universal_KEGG_dict)
            EuPath_reaction_gene_dict[hf.transcript_to_gene_id(gene_id)] = [r for r in RXNs]
            add_EuPath.append(RXNs)
file1.close()

# remove genes with no reaction mapping
reaction_gene_dict_keep = dict()
for gene, rxn_list in EuPath_reaction_gene_dict.items():
    if len(rxn_list) > 0: reaction_gene_dict_keep[gene] = rxn_list

# flip dictionary to focus on reactions rather than genes
new_dict = dict()
for gene, rxn_list in reaction_gene_dict_keep.items():
    for rxn in rxn_list:
        if rxn in new_dict.keys():
            gene_list = new_dict[rxn]
            new_dict[rxn].extend([gene]) #([gene+'-t36_1-p1'])
        else:
            new_dict[rxn] = [gene] #[gene+'-t36_1-p1']
add_to_model_from_EuPathDB = new_dict

# need to move reactions in add_to_model to appropriate compartments
# add reactions with gene association]
in_model_already_reactions = list()
in_model_already_genes = list()

for rxn_id, gene_list in add_to_model_from_EuPathDB.items():
    
    # prep genes
    if rxn_id in [r.id for r in new_model.reactions]: 
        in_model_already_reactions.append(rxn_id)
        rxn_not_present_yet = False
    else: rxn_not_present_yet = True
        
    for gene_id in gene_list: 
        if sum([gene_id in x for x in [g.id for g in new_model.genes]]) > 0: #gene_id in [g.id for g in new_model.genes]: 
            in_model_already_genes.append(gene_id)
            
    # add reaction if not in model yet
    if rxn_not_present_yet:
        if len(gene_list) > 1: use_gene = ' or '.join(gene_list)
        else: use_gene = gene_list[0]
        rxn = universal_model.reactions.get_by_id(rxn_id)
        for met in rxn.metabolites:
            if met.id not in new_model.metabolites:
                new_model.add_metabolites([met])
        rxn.notes['EuPathDB'] = 'version 50'
        new_model.add_reactions([rxn])
        new_model.reactions.get_by_id(rxn_id).gene_reaction_rule = use_gene
    # add gene to reaction if rxn is already in model yet
    else:
        existing_genes = [g.id for g in new_model.reactions.get_by_id(rxn_id).genes]
        for gene_id in gene_list:
            if sum([gene_id in x for x in [g.id for g in new_model.genes]]) == 0:
                existing_genes.extend([gene_id])
        unique_gene_list = list(set(existing_genes))
        new_model.reactions.get_by_id(rxn_id).gene_reaction_rule = ' or '.join(unique_gene_list)
        
# save genes that were on EuPathDB but not detected with Diamond
EuPath_only = list()
for key in reaction_gene_dict_keep.keys():
    if list(filter(lambda x: x.startswith(key), [gene.id for gene in g])):
        EuPath_only.append(key)
logging.info('genes added from EuPathDB:')
logging.info(EuPath_only)

row_index =0 #= SPECIES_ID
modifications.loc[row_index,'genes_only_on_EuPathDB'] = len(EuPath_only)


# moving reactions based on compartment
model = new_model
logging.info('finding good or bad reactions')
good_rxns, bad_rxns = hf.id_bad_compartment_rxns(model,compartment,compartment_options)
# get reactions that use/make at least one metabolite that is in an inappropariate compartment

logging.info('found good or bad reactions, now doing things')
remove_rxn, add_reaction, bad_rxns_keep_rewrite = hf.move_bad_rxns(model,bad_rxns,universal_dict_with_alts, compartment)
if (len((bad_rxns)) == (len((remove_rxn)) + len((bad_rxns_keep_rewrite)))):
    logging.info('bad reactions are split into remove reactions and bad reactions to rewrite - math is good')

logging.info('no. reactions removed')
logging.info(len(remove_rxn))
logging.info('no. reactions to add')
logging.info(len(add_reaction))

# remove reactions
x = len(model.reactions)
y = len(model.metabolites)
model.remove_reactions(remove_rxn)
model.repair()

if len(remove_rxn) != (x - len(model.reactions)):
    logging.info('error - reaction not removed properly')
x1 = x - len(model.reactions)
y1 = y - len(model.metabolites)
x1_2 = len(model.reactions)

# save this number
inappropriate_compartments_that_remain = (len(bad_rxns_keep_rewrite)/len(model.reactions))*100

row_index =0 #= SPECIES_ID
modifications.loc[row_index,'species'] = SPECIES_ID
modifications.loc[row_index,'reactions_removed1'] = x1

for rxn_id in add_reaction: #there are ids in add_reaction that are in the model already
    rxn = universal_model.reactions.get_by_id(rxn_id).copy()
    for met in rxn.metabolites:
        if met.id not in [m.id for m in model.metabolites]:
            model.add_metabolites(met.copy())
rxns_to_add_list = [universal_model.reactions.get_by_id(x).copy() for x in add_reaction if x not in [r.id for r in model.reactions]]
# if reaction is already there, it is because the reaction was in multiple compartments
model.add_reactions(rxns_to_add_list)
row_index = 0#modifications.species == SPECIES_ID
modifications.loc[row_index,'reactions_added'] = len(model.reactions) - x1_2

# make sure all reactions can carry flux, problem with some versions of the universal model
for rxn in model.reactions:
    if rxn.lower_bound == 0 and rxn.upper_bound == 0:
        logging.info(rxn.id + ' has bounds == 0 in '+key)
        rxn.lower_bound = -1000.
        rxn.upper_bound = 1000.
        # NOTHING SHOULD PRINT - this was a problem in CarveMe

fix_these_reactions = list(set([model.reactions.get_by_id(x) for x in bad_rxns_keep_rewrite]))

og = len(model.reactions)
og_mets= len(model.metabolites)

logging.info('starting to move reactions to the right compartment')
model, error_dict_to_print, reactions_added, transport_for_inappropariate_compartment = hf.fixing_reaction_compartment(fix_these_reactions,model)
logging.info(error_dict_to_print) 

model, unused = hf.prune_unused_metabolites2(model)
logging.info('finished moving reactions to the right compartment')

# ADD IN THIS WARNING:
for c in list(set(compartment_options) - set(compartment)):
    if c in ([hf.get_comp(model,x.id) for x in model.metabolites]):
        logging.info('ERROR, UNACCEPTABLE COMPARTMENTS:')
        logging.info(c)

logging.info('reactions added overall:')
logging.info(len(model.reactions) - og)
transport_for_inappropariate_compartment_dict = list(set(transport_for_inappropariate_compartment))

l2 = list()
for rxn in model.reactions:
    for suffix in [hf.get_comp(model,m.id) for m in rxn.metabolites]:
        l2.append(suffix)
logging.info('compartments:')
logging.info(set(l2))

os.chdir(data_path)
modifications.to_csv('model_modifications/model_modifications_'+SPECIES_ID+day+'.csv')
with open("./percent_wrong_comp/percent_reactions_in_wrong_compartment_"+SPECIES_ID+day+".csv", "w") as text_file:
    text_file.write(str(inappropriate_compartments_that_remain))

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
cobra.io.save_json_model(model, "final_denovo_"+SPECIES_ID+".json")
