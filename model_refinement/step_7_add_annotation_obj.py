import cobra
import re
import os
import glob
from cobra import Model, Reaction, Metabolite
import logging
from datetime import datetime

log_file_path = "/home/mac9jc/paradigm/model_generation_logs"
model_path = "/home/mac9jc/paradigm/models"

os.chdir(log_file_path)
day = datetime.now().strftime('%d_%m_%Y')
logging.basicConfig(filename='add_annotation.log'.format(day), level=logging.INFO, filemode='w')
logger = logging.getLogger(__name__)

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

compartment_dict = {'c': 'cytoplasm', 'e': 'extracellular', 'm': 'mitochondrion', 'fv': 'food vacuole', 'ap':'apicoplast','k':'kinetoplast','glc':'glycosome'}
version = '44'

os.chdir(model_path)
for file in glob.glob("gf_*.xml"):

    SPECIES_ID = file.split('gf_')[1] # ID is annotation filename minus directory
    SPECIES_ID = SPECIES_ID.split('.json')[0] # get rid of extension
    if SPECIES_ID.startswith('no_ortho'):
        SPECIES_ID = SPECIES_ID.split('no_ortho_')[1] # get rid of no_ortho

    if SPECIES_ID in plasmodb: # Plasmodium = cytosol, extracellular, mitochondrdia, apicoplast, food vacuole
        compartment = ["c","e","m","ap","fv"]
    elif SPECIES_ID in tritrypdb: # Leishmania = cytosol, extracellular, mitochondrdia, kinetoplast, glycosome
        compartment = ["c","e","m","k","glc"]
    elif SPECIES_ID in cryptodb: # Cryptosporidium = cytosol, extracellular, pseudomitochondria (USE MITO)
        compartment = ["c","e","m"]
    elif SPECIES_ID in toxodb: # Toxoplasma = cytosol, extracellular, mitochondrdia, apicoplast
        compartment = ["c","e","ap","m"]
    else: compartment = ["c","e"]
    
    model = cobra.io.read_sbml_model(file)

    # add exchange for all extracellular mets
    for met in model.metabolites:
        if met.id.endswith('_e'):
            if 'EX_'+met.id not in [r.id for r in model.reactions]:
                model.add_boundary(met, type="exchange", lb = -1000., ub = 1000.)

    compartment_use = {key:value for key, value in compartment_dict.items() if key in compartment}
    model.compartments = compartment_use
    model.repair()

    logger.info(file)    

    if 'no_ortho_P' in file:
        model.notes = 'This reconstruction represents the metabolism of {} and was generated as part of ParaDIGM, v1. it was NOT semi-curated using orthology.'.format(SPECIES_ID)
    elif SPECIES_ID in plasmodb:
        model.notes  = 'This reconstruction represents the metabolism of {} and was generated as part of ParaDIGM, v1. it was semi-curated using orthology.'.format(SPECIES_ID)
    else: model.notes = 'This reconstruction represents the metabolism of {} and was generated as part of ParaDIGM, v1.'.format(SPECIES_ID)
    model.annotation["taxonomy"] = "must add this ID"

    if SPECIES_ID in plasmodb: database = 'PlasmoDB'
    elif SPECIES_ID in tritrypdb: database = 'TriTrypDB'
    elif SPECIES_ID in cryptodb: database = 'CryptoDB'
    elif SPECIES_ID in toxodb: database = 'ToxoDB'
    elif SPECIES_ID in giardiadb: database = 'GiardiaDB'
    elif SPECIES_ID in piroplasmadb: database = 'PiroplasmaDB'
    elif SPECIES_ID in microsporidiadb: database = 'MicrosporidiaDB'
    elif SPECIES_ID in trichdb: database = 'TrichDB'
    elif SPECIES_ID in amoebadb: database = 'AmoebaDB'
    else: database = 'ERROR'
    model.annotation["genome"] = "https://{}.org/common/downloads/release-{}/{}/fasta/data/{}-{}_{}_Genome.fasta".format(database.lower(),version,SPECIES_ID,database,version,SPECIES_ID)
    logger.info(model.annotation["genome"])

    model.annotation["DOI"] = "pending"
    model.annotation["authors"] = ["Maureen Carey"]
    model.annotation["corresponding_author"] = "Maureen Carey, mac9jc@virginia.edu"

    temp_species = re.sub( r"([A-Z])", r" \1", SPECIES_ID).split()[0]
    model.annotation["species"] = temp_species[0]+'. '+temp_species[1:]   
    #model.annotation["species"] = re.sub( r"([A-Z])", r" \1", SPECIES_ID).split()[0]
    model.annotation["strain"] = ''.join(re.sub( r"([A-Z])", r" \1", SPECIES_ID).split()[1:])
    model.annotation["tissue"] = "NA"
    logger.info(model.annotation["species"])
    logger.info(model.annotation["strain"])

    model.id = 'i{}_v1'.format(temp_species[0:4])
    model.name = 'i{}'.format(temp_species[0:4])
    logger.info(model.name)

    model.annotation["terms_of_distribution"] = "CC-BY"
    model.annotation["updated"] = day
    
    cobra.io.write_sbml_model(model, "final_{}.xml".format(SPECIES_ID))
    logger.info('done with {}'.format(file))
