import cobra
import pandas as pd
import os
import subprocess
import glob
import json
from cobra import Model, Reaction, Metabolite
import helper_functions_1 as hf
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


day = datetime.now().strftime('%d_%m_%Y')
logging.basicConfig(filename='step3_{}_{}.log'.format(SPECIES_ID,day), level=logging.INFO, filemode='w')
logger = logging.getLogger(__name__)

logging.info('BEGIN STEP 3')
logging.info(SPECIES_ID)

data_path = "/home/mac9jc/paradigm/data"
model_path = "/home/mac9jc/paradigm/models"
og_path = "/home/mac9jc/paradigm/"

# universal reaction bag for model generation
os.chdir(model_path)
universal_model = cobra.io.load_json_model('universal_model_oct26_2018.json')
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

# get metabolites involved in each reaction in universal model in dictionary format
universal_dict = dict() 
for rxn in universal_model.reactions:
    rxn_dict = dict()

    check_rxn_products = rxn.products
    check_rxn_reactants = rxn.reactants
    all_mets = rxn.metabolites
    compart = set([x.id[-2:] for x in all_mets])

    rxn_dict['reactants'] = [x.id[:-2] for x in check_rxn_reactants]
    rxn_dict['products'] = [x.id[:-2] for x in check_rxn_products]
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

# database mapping - this is for release 39, must update if using a different EuPathDB release
plasmodb = ["PadleriG01","PbergheiANKA","PbillcollinsiG01","PblacklockiG01","Pchabaudichabaudi","PcoatneyiHackeri","PcynomolgiB","PcynomolgiM","Pfalciparum3D7","Pfalciparum7G8","PfalciparumCD01","PfalciparumDd2","PfalciparumGA01","PfalciparumGB4","PfalciparumGN01","PfalciparumHB3","PfalciparumIT","PfalciparumKE01","PfalciparumKH01","PfalciparumKH02","PfalciparumML01","PfalciparumSD01","PfalciparumSN01","PfalciparumTG01","PfragileNilgiri","PgaboniG01","PgaboniSY75","Pgallinaceum8A","PinuiSanAntonio1","PknowlesiH","PknowlesiMalayanPk1A","PmalariaeUG01","PovalecurtisiGH01","PpraefalciparumG01","PreichenowiCDC","PreichenowiG01","PrelictumSGS1-like","PvinckeipetteriCR","Pvinckeivinckeivinckei","PvivaxP01","PvivaxSal1","Pyoeliiyoelii17X","Pyoeliiyoelii17XNL","PyoeliiyoeliiYM"]
cryptodb = ["Candersoni30847","Chominis30976","ChominisTU502","ChominisTU502_2012","ChominisUdeA01","CmeleagridisUKMEL1","CmurisRN66","CparvumIowaII","CtyzzeriUGA55", "Cubiquitum39726","CveliaCCMP2878", "GniphandrodesUnknown", "VbrassicaformisCCMP3155"]
giardiadb = ["GintestinalisAssemblageADH", "GintestinalisAssemblageAWB", "GintestinalisAssemblageBGS", "GintestinalisAssemblageBGS_B", "GintestinalisAssemblageEP15", "SsalmonicidaATCC50377"]
tritrypdb = ["BayalaiB08-376","CfasciculataCfCl","EmonterogeiiLV88","LaethiopicaL147", "LarabicaLEM1108", "LbraziliensisMHOMBR75M2903", "LbraziliensisMHOMBR75M2904", "LdonovaniBPK282A1", "LenriettiiLEM3045", "LgerbilliLEM452","LinfantumJPCM5", "LmajorFriedlin", "LmajorLV39c5", "LmajorSD75.1", "LmexicanaMHOMGT2001U1103", "LpanamensisMHOMCOL81L13","LpanamensisMHOMPA94PSC1", "LpyrrhocorisH10", "LseymouriATCC30220", "LspMARLEM2494", "LtarentolaeParrotTarII", "LtropicaL590", "LturanicaLEM423", "PconfusumCUL13","TbruceigambienseDAL972", "TbruceiLister427", "TbruceiTREU927", "TcongolenseIL3000", "TcruziCLBrener", "TcruziCLBrenerEsmeraldo-like", "TcruziCLBrenerNon-Esmeraldo-like", "TcruzicruziDm28c","TcruziDm28c", "TcruzimarinkelleiB7", "TcruziSylvioX10-1", "TcruziSylvioX10-1-2012","TevansiSTIB805", "TgrayiANR4", "TrangeliSC58", "TvivaxY486", "TtheileriEdinburgh"]
# MUST ASK PERMISSION TO USE MANY OF THE TRITRYP GENOMES - S.M. Beverley at Wash U  as of Oct 4 2018
trichdb = ["TvaginalisG3"]
amoebadb = ["AcastellaniiNeff", "EdisparSAW760", "EhistolyticaHM1IMSS-A", "EhistolyticaHM1IMSS-B", "EhistolyticaHM1IMSS", "EhistolyticaHM3IMSS", "EhistolyticaKU27", "EinvadensIP1", "EmoshkovskiiLaredo", "EnuttalliP19", "NfowleriATCC30863"]
toxodb = ["CcayetanensisCHN_HEN01", "CsuisWienI","EacervulinaHoughton", "EbrunettiHoughton", "EfalciformisBayerHaberkorn1970", "EmaximaWeybridge", "EmitisHoughton", "EnecatrixHoughton", "EpraecoxHoughton", "EtenellaHoughton", "HhammondiHH34", "NcaninumLIV", "SneuronaSN3", "SneuronaSOSN1", "TgondiiARI", "TgondiiFOU", "TgondiiGAB2-2007-GAL-DOM2", "TgondiiGT1", "TgondiiMAS", "TgondiiME49", "Tgondiip89", "TgondiiRH", "TgondiiRUB", "TgondiiTgCatPRC2", "TgondiiVAND", "TgondiiVEG"]
microsporidiadb = ["AalgeraePRA109", "AalgeraePRA339", "AspWSBS2006","EaedisUSNM41457", "EbieneusiH348", "EcanceriGB1","EcuniculiEC1", "EcuniculiEC2", "EcuniculiEC3","EcuniculiEcunIII-L","EcuniculiGBM1", "EhellemATCC50504", "EhellemSwiss", "EhepatopenaeiTH1","EintestinalisATCC50506", "EromaleaeSJ2008","Heriocheircanceri","HeriocheirGB1", "MdaphniaeUGP3", "NausubeliERTm2", "NausubeliERTm6", "NbombycisCQ1", "NceranaeBRL01", "NceranaePA08_1199","NdisplodereJUm2807","NparisiiERTm1", "NparisiiERTm3", "OcolligataOC4", "PneurophiliaMK1", "Slophii42_110", "ThominisUnknown", "VcorneaeATCC50505", "Vculicisfloridensis"]
piroplasmadb = ["BbigeminaBOND", "BbovisT2Bo", "BmicrotiRI","BovataMiyake", "CfelisWinnie", "TannulataAnkara", "TequiWA", "TorientalisShintoku", "TparvaMuguga"]
fungidb = ["AaculeatusATCC16872", "AbrasiliensisCBS101740", "AcampestrisIBT28561", "Acandida2VRR", "AcarbonariusITEM5010", "AclavatusNRRL1", "AfischeriNRRL181","AflavusNRRL3357","AfumigatusA1163","AfumigatusAf293", "AglaucusCBS516.65","AinvadansNJM9701", "AlaibachiiNc14", "AluchuensisCBS106.47", "AmacrogynusATCC38327", "AnidulansFGSCA4", "AnigerATCC1015", "AnigerCBS513-88", "AnovofumigatusIBT16806", "AochraceoroseusIBT24754","AoryzaeRIB40", "AsteyniiIBT23096", "AsydowiiCBS593.65", "AterreusNIH2624", "AtubingensisCBS134.48", "AversicolorCBS583.65", "AwentiiDTO134E9","AzonataCBS506.65", "BcinereaB05-10", "BdendrobatidisJEL423", "CalbicansSC5314", "CalbicansSC5314_B", "CalbicansWO1", "CaurisB8441", "Ccinereaokay7-130", "CdeuterogattiiR265", "CgattiiCA1873", "CgattiiEJB2", "CgattiiIND107", "CgattiiWM276", "CglabrataCBS138", "CimmitisH538-4", "CimmitisRS", "ClusitaniaeATCC42720", "CneoformansB-3501A", "CneoformansH99", "CneoformansJEC21", "CneoformansKN99", "CposadasiiC735deltSOWgp", "CposadasiiRMSCC3488", "CposadasiiRMSCC3700", "CposadasiiSilveira", "FfujikuroiIMI58289", "FgraminearumPH-1", "Foxysporum26406", "Foxysporum4287", "Foxysporum54006", "FoxysporumFo47", "Foxysporumrace1", "Foxysporumrace4", "Fverticillioides7600", "HarabidopsidisEmoy2", "HcapsulatumG186AR", "HcapsulatumG217B", "HcapsulatumH143", "HcapsulatumH88", "HcapsulatumNAm1", "McircinelloidesCBS277-49", "MglobosaCBS7966", "Mlarici-populina98AG31", "Moryzae70-15", "MoryzaeBR32", "NcrassaOR74A", "NdiscretaFGSC8579", "NtetraspermaFGSC2508", "PaphanidermatumDAOMBR444", "ParrhenomanesATCC12531", "PblakesleeanusNRRL1555", "PbrasiliensisPb03", "PbrasiliensisPb18", "PcapsiciLT1534", "PchrysosporiumRP-78", "PcinnamomiCBS144-22", "PgraminisCRL75-36-700-3", "PinfestansT30-4", "PirregulareDAOMBR486", "PiwayamaiDAOMBR242034", "PjiroveciiSE8", "PlutziiPb01", "PparasiticaINRA-310", "PramorumPr-102", "PrubensWisconsin54-1255", "PsojaeP6497", "PultimumBR650", "PultimumDAOMBR144", "PvexansDAOMBR484", "RdelemarRA99-880", "ScerevisiaeS288c", "SdiclinaVS20", "SjaponicusyFS275", "Smacrosporak-hell", "SoctosporusyFS286", "SparasiticaCBS223", "Spombe972h", "SpunctatusDAOMBR117", "SreilianumSRZ2", "Sschenckii1099-18", "Ssclerotiorum1980UF-70", "TmarneffeiATCC18224", "TmesentericaDSM1558", "TreeseiQM6a", "TstipitatusATCC10500", "Umaydis521", "Ureesii1704", "YlipolyticaCLIB122", "ZtriticiIPO323","AkawachiiIFO4308", "AnigerN402ATCC64974","CgattiiNT10","CparapsilosisCDC317","FproliferatumET1", "LprolificansJHH5317", "SapiospermumIHEM14462","Sbrasiliensis5110","YlipolyticaCLIB89W29"]

# load annotation data file
annotations_dict = dict()
columns = ['query_gene', 'BiGG_gene', 'pident', 'length', 'mismatch', 'gapopen','qstart', 'qend', 'sstart', 'send', 'evalue', 'score']

# modified for Rivanna: read in the annotation tsv
annotations_dict = {}
annotations_dict[SPECIES_ID] = pd.read_csv(annotation_fname,'\t')
columns = ['query_gene', 'BiGG_gene', 'pident', 'length', 'mismatch', 'gapopen','qstart', 'qend', 'sstart', 'send', 'evalue', 'score']
annotations_dict[SPECIES_ID].columns = columns

logging.info('loading gprs')
# map annotations to BiGG functions for model generation
os.chdir(data_path)
gprs = pd.read_csv('./bigg_gprs.csv') # From CarveMe
gprs.reaction = [x[2:] for x in gprs.reaction]
gprs = gprs[gprs.reaction.isin([rxn.id for rxn in universal_model.reactions])] # updated from CarveMe

logging.info('loaded gprs')

def merge_subunits(genes): # From CarveMe
    """ Merge list of protein subunit genes into complex
    Args: genes (pandas.Series): list of genes
    Returns: str: boolean rule
    """
    genes = genes.dropna()

    if len(genes) == 0:
        return None
    else:
        protein = ' and '.join(sorted(genes))
        if len(genes) > 1:
            return '(' + protein + ')'
        else:
            return protein

def merge_subunit_scores(scores): # From CarveMe
    """ Merge scores of all genes in a protein complex.
    Calculates the mean score among all subunits.
    Args: scores: individual gene scores
    Returns: float: merged score
    """
    return scores.fillna(0).mean()

def merge_proteins(proteins): # From CarveMe
    """ Merge all isozymes that catalyze a given reaction.
    Automatically removes all isozymes with missing score.
    Args: proteins (pandas.Series): list of proteins
    Returns: str: boolean rule
    """
    proteins = set(proteins.dropna())
    if not proteins:
        return None
    gpr_str = ' or '.join(sorted(proteins))
    if len(proteins) > 1:
        return '(' + gpr_str + ')'
    else:
        return gpr_str

def merge_protein_scores(scores): # From CarveMe
    """ Merge scores of all isozymes that catalyze a given reaction.
    Calculates the maximum score among all isozymes.
    Args: scores (pandas.Series): protein scores
    Returns: float: merged score
    """
    return scores.max(skipna=True)

def reaction_scoring(annotation, gprs, spontaneous_score=0.0, debug_output=None): # From CarveMe
    """ Calculate reaction scores using new eggnog output.
    Args: annotation (pandas.DataFrame): gene annotation results
        gprs (pandas.DataFrame): BiGG GPR rules
        spontaneous_score (float): score to give to spontaneous reactions (default: 0.0)
    Returns: pandas.DataFrame: reaction scores
    """

    # filter best match for each gene
    gene2gene = annotation.sort_values(by='score', ascending=False) \
                          .groupby('BiGG_gene', as_index=False).apply(lambda x: x.iloc[0])
    # merge with gpr table
    gprs['BiGG_gene'] = gprs.apply(lambda row: '{}.{}'.format(row['model'], row['gene'][2:]), axis=1)
    gene_scores = pd.merge(gene2gene, gprs, how='right')
    # add default scores for spontaneous genes
    spontaneous = {'G_s0001', 'G_S0001', 'G_s_0001', 'G_S_0001', 'G_KPN_SPONT'}
    gene_scores.loc[gene_scores.gene.isin(spontaneous), 'score'] = spontaneous_score
    gene_scores.loc[gene_scores.gene.isin(spontaneous), 'query_gene'] = 'spontaneous'
    # from gene to protein scores
    protein_scores = gene_scores.groupby(['protein', 'reaction', 'model'], as_index=False) \
        .agg({'query_gene': merge_subunits, 'score': merge_subunit_scores})
    protein_scores.rename(columns={'query_gene': 'GPR'}, inplace=True)
    # from protein to reaction scores
    reaction_scores = protein_scores.groupby(['reaction'], as_index=False) \
        .agg({'GPR': merge_proteins, 'score': merge_protein_scores}).dropna()
    return(reaction_scores)

scores_dict2 = dict()
for species, annotations in annotations_dict.items():
    reaction_scores = reaction_scoring(annotations, gprs)
    # scores = dict(reaction_scores[['reaction', 'normalized_score']].values)
    scores_dict2[species] = reaction_scores #scores
    # carveme will maximize positive scores and minimize negative scores while maintaining a functional network
logging.info("done with scoring")

# make model from universal model
new_model_dict = dict()
for x in scores_dict2.keys():
    logging.info(x)
    
    if '.tsv' in x:
        logging.info('.tsv in the annotations file string, this might cause problems')
        species = x.split('.tsv')[0]
    else:
        species = x
    if '_BiGG' in x:
        logging.info('_BiGG in the annotations file string, this might cause problems')
        species = species.split('_BiGG')[0]

    new_model_dict[species] = universal_model.copy()
    logging.info('copied universal')
    new_model_dict[species].name = species
    new_model_dict[species].id = species

    logging.info(species)
    starting = len(new_model_dict[species].reactions)

    keep_scores = scores_dict2[x].loc[scores_dict2[x].score>10]
    if len(keep_scores.reaction) == len(set(keep_scores.reaction)):

        rxns_to_add = dict()
        scores_for_rxns = dict()
        for index, row in keep_scores.iterrows():
            rxns_to_add[row['reaction']] = row['GPR']
            scores_for_rxns[row['reaction']] = row['score']

        new_model_dict[species].remove_reactions([rxn for rxn in new_model_dict[species].reactions if rxn.id not in rxns_to_add.keys()])

        if not [rxn.id for rxn in new_model_dict[species].reactions if rxn.gene_reaction_rule != '']:
            for rxn in new_model_dict[species].reactions:
                    if rxn.gene_reaction_rule == '':
                        new_model_dict[species].reactions.get_by_id(rxn.id).gene_reaction_rule = rxns_to_add[rxn.id]
                        new_model_dict[species].reactions.get_by_id(rxn.id).notes['CarveMe score'] = {rxns_to_add[rxn.id] : scores_dict2[x].loc[scores_dict2[x].GPR == rxns_to_add[rxn.id]]['score']}
        else:
            logging.info('some reactions already have GPRs')

        logging.info('made new model from universal')
        new_model_dict[species].repair()

        if len(rxns_to_add.keys()) == len(new_model_dict[species].reactions):
            if starting > len(rxns_to_add.keys()):
                logging.info(' ')
            else:
                logging.info('error with original model, reactions already removed')
        else:
            logging.info('error with reaction removal, resultant len(model.reactions) != rxns_to_keep')
    else:
        logging.info('duplicate keep_scores.reaction')

    if len(rxns_to_add.keys()) != len(new_model_dict[species].reactions):
        logging.info('error in universal reaction pruning')
    model2 = new_model_dict[species].copy()
    logging.info('made DIY1')
    logging.info(len(model2.reactions))
    os.chdir(model_path)
    cobra.io.save_json_model( model2, "DIY1_"+species+".json")
    logging.info('saved DIY1')
#temp
    logging.info('notes field')
    logging.info([rxn.notes for rxn in model2.reactions])
    logging.info('------------------------------------------')
    
# remove duplciate reactions in mulitple compartments
total_compartments = ["_c","_e","_m","_ap","_fv","_k","_glc","_pm"]
# cytosol, extracellular, mitochondrdia, apicoplast, food vacuole, kinetoplast, glycosome, pseudomitochondria

compartment_dictionary = dict()
for species in new_model_dict.keys():

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

columns = ['species','reactions_removed1','mets_removed1','reactions_removed2','mets_removed2','reactions_added', 'mets_added','gene_change']
modifications = pd.DataFrame(index = new_model_dict.keys(), columns=columns)
inappropriate_compartments_that_remain = dict()
transport_for_inappropariate_compartment_dict = dict()

for species, model in new_model_dict.items():

    logging.info(species)
    logging.info('finding good or bad reactions')
    compartment = compartment_dictionary[species]
    not_compartments = compartment_options - set(compartment)

    # get reactions that use/make at least one metabolite that is in an inappropariate compartment
    good_rxns = list()
    bad_rxns = list()
    not_compartments = [x+' ' for x in not_compartments]
    for rxn_object in model.reactions: # if a reaction does not contain any bad compartments
        rxn_bad_counter = 0
        for x in not_compartments:
            if x in rxn_object.reaction or rxn_object.reaction.endswith(x[:-1]):
                rxn_bad_counter = rxn_bad_counter + 1
        if rxn_bad_counter == 0:
            good_rxns.append(rxn_object.id)
        else:
            bad_rxns.append(rxn_object.id)

    logging.info('found good or bad reactions, now doing things')
    bad_rxns_keep_rewrite = list()
    add_reaction = list()
    remove_rxn = list()
    for rxn_id in bad_rxns:
        if len(universal_dict_with_alts[rxn_id]['alternative_reactions']) == 0:
            bad_rxns_keep_rewrite.append(rxn_id) # no alternative, keep reaction - will have to change via strings
        else:
            alt_rxns = universal_dict_with_alts[rxn_id]['alternative_reactions']
            keep_og = 0
            for alt_rxn_1, locations in alt_rxns.items():
                keep_alt = 0
                for loc in locations:
                    if loc in compartment: keep_alt = keep_alt
                    else: keep_alt = keep_alt + 1
                if keep_alt == 0:
                    keep_og = 1
                    add_reaction.append(alt_rxn_1)
                else:
                    keep_og = keep_og
            if keep_og == 0:
                bad_rxns_keep_rewrite.append(rxn_id) # no usable alternative - will have to change via strings
            else:
                remove_rxn.append(rxn_id)
    add_reaction = list(set(add_reaction))

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
        logging.info('reaction not removed properly')
    x1 = x - len(model.reactions)
    y1 = y - len(model.metabolites)

    # save this number
    inappropriate_compartments_that_remain[species] = (len(bad_rxns_keep_rewrite)/len(model.reactions))*100
    logging.info((len(bad_rxns_keep_rewrite)/len(model.reactions))*100)

    modifications.species.loc[species] = species
    modifications.reactions_removed1.loc[species] = x1
    #     modifications.mets_removed1.loc[species] = y1

    for rxn_id in add_reaction: #there are ids in add_reaction that are in the model already
        rxn = universal_model.reactions.get_by_id(rxn_id).copy()
        for met in rxn.metabolites:
            if met.id not in [m.id for m in model.metabolites]:
                model.add_metabolites(met.copy())
    rxns_to_add_list = [universal_model.reactions.get_by_id(x).copy() for x in add_reaction if x not in [r.id for r in model.reactions]]
    # if reaction is already there, it is because the reaction was in multiple compartments
    # logging.info([rxn.id for rxn in rxns_to_add_list if rxn.id not in add_reaction])
    model.add_reactions(rxns_to_add_list)
    modifications.reactions_added.loc[species] = len(model.reactions) - x1 # CHECK

#     if '_x' in ([x.id[-2:] for x in model.metabolites]): logging.info('ERROR, UNACCEPTABLE COMPARTMENTS, should be some here')
#     logging.info('the right number of reactions are being added and removed')

    for rxn in model.reactions:
        if rxn.lower_bound == 0 and rxn.upper_bound == 0:
            logging.info(rxn.id + ' has bounds == 0 in '+key)
            rxn.lower_bound = -1000.
            rxn.upper_bound = 1000.
            # NOTHING SHOULD PRINT - this was a problem in CarveMe

    new_model_dict[species] = model
    os.chdir(model_path)
    cobra.io.save_json_model(model, "DIY2_"+species+".json")
    ### UPDATE NAMING HERE

    fix_these_reactions_list = list(set([model.reactions.get_by_id(x) for x in bad_rxns_keep_rewrite]))

    reactions_added = list()
    transport_for_inappropariate_compartment = list()

    og = len(model.reactions)
    og_mets= len(model.metabolites)

    logging.info('starting to move reactions to the right compartment, if this function isnt yet in BiGG')
    for rxn in fix_these_reactions_list:

        if [hf.met_ids_without_comp(model,x.id) for x in rxn.reactants] == [hf.met_ids_without_comp(model,x.id) for x in rxn.products]:
            # remove things like x_p + y_p => x_e + y_e
            transport_for_inappropariate_compartment.append(rxn.id)
            new_rxn = list()
        else:

            new_rxn = Reaction()
            met_dict = dict()

            for met in rxn.metabolites:

                if hf.get_comp(model,met.id) == '_p': # move periplasmic metabolites to extracellular instead of cytosol
                    if hf.met_ids_without_comp(model,met.id)+'_e' not in [x.id for x in model.metabolites]:
                        met2 = met.copy()
                        met2.id = hf.met_ids_without_comp(model,met.id)+'_e'
                        met_dict[met2] = rxn.metabolites[met]
                        model.add_metabolites(met2) # []
                    else:
                        met2 = model.metabolites.get_by_id(hf.met_ids_without_comp(model,met.id)+'_e')
                        met_dict[met2] = rxn.metabolites[met]
                else: # non periplasmic metabolite
                    if hf.met_ids_without_comp(model,met.id)+'_c' not in [x.id for x in model.metabolites]:
                        met2 = met.copy()
                        met2.id = hf.met_ids_without_comp(model,met.id)+'_c'
                        met_dict[met2] = rxn.metabolites[met]
                        model.add_metabolites(met2) # []
                    else:
                        met2 = model.metabolites.get_by_id(hf.met_ids_without_comp(model,met.id)+'_c')
                        met_dict[met2] = rxn.metabolites[met]

        # fix reaction variables
        if new_rxn:
            new_rxn.add_metabolites(met_dict)
            new_rxn.name = rxn.name
            new_rxn.id = rxn.id+'c'
            new_rxn.lower_bound = rxn.lower_bound
            new_rxn.upper_bound = rxn.upper_bound
            new_rxn.gene_reaction_rule = rxn.gene_reaction_rule
            model.add_reactions([new_rxn])
            reactions_added.append(new_rxn.id)

        model.remove_reactions([rxn])

    model.repair()
    logging.info('finished moving reactions to the right compartment')

    logging.info('reactions added overall')
    logging.info(len(model.reactions) - og)
    # logging.info('mets added')
    # logging.info(len(model.metabolites) - og_mets)
    transport_for_inappropariate_compartment_dict[species] = list(set(transport_for_inappropariate_compartment))

    new_model_dict[species] = model

    os.chdir(model_path)
    cobra.io.save_json_model(model, "DIY3_"+species+".json")
    l2 = list()
    for rxn in model.reactions:
        for suffix in [m.id[-2:] for m in rxn.metabolites]:
            l2.append(suffix)
    logging.info('compartments:')
    logging.info(set(l2))

os.chdir(data_path)
modifications.to_csv('model_modifications_'+SPECIES_ID+'Dec.csv')
pd.DataFrame.from_dict(inappropriate_compartments_that_remain, orient="index").to_csv("./percent_reactions_in_wrong_compartment_"+SPECIES_ID+"_Dec.csv")

def prune_unused_metabolites2(cobra_model):
    """ USE THIS UNTIL AUG 31 UPDATES ARE INTEGRATED INTO MASTER COBRAPY BRANCH
    Remove metabolites that are not involved in any reactions and
    returns pruned model
    Parameters
    ----------
    cobra_model: class:`~cobra.core.Model.Model` object
        the model to remove unused metabolites from
    Returns
    -------
    output_model: class:`~cobra.core.Model.Model` object
        input model with unused metabolites removed
    inactive_metabolites: list of class:`~cobra.core.reaction.Reaction`
        list of metabolites that were removed
    """
    output_model = cobra_model.copy()
    inactive_metabolites = [m for m in output_model.metabolites if len(m.reactions) == 0]
    output_model.remove_metabolites(inactive_metabolites)
    return output_model, inactive_metabolites

for species, model in new_model_dict.items():
    
    model_pruned, unused = prune_unused_metabolites2(model)
    model_pruned.solver = 'gurobi'
    new_model_dict[species] = model_pruned
    
    # check compartments
    list_om= list()
    list_om2= list()
    for rxn in model_pruned.reactions:
        for m in rxn.metabolites:
            list_om.append(m.id[-2:])
    for m in model_pruned.metabolites:
        list_om2.append(m.id[-2:])
    if set(list_om) != set(list_om2):
        logging.info('extra compartments are present, pruning of unused metabolites did not work')
    
    os.chdir(model_path)
    cobra.io.save_json_model(model_pruned, "final_denovo_"+species+".json")
