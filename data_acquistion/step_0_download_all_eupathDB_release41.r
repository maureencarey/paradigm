# library(tidyverse)
library(RCurl)
# Feb 9th 2018 MAC, updated May 23rd 2018, MAC, updated July 24th 2018 MAC, updated for release 39 Oct 4th, 2018 MAC, updated for release 41 on Dec 11, 2018

version = 41
database_lc = c("plasmodb","cryptodb","giardiadb","tritrypdb","trichdb","amoebadb","toxodb","microsporidiadb", "piroplasmadb")
DataBase_uc = c("PlasmoDB","CryptoDB","GiardiaDB","TriTrypDB","TrichDB","AmoebaDB","ToxoDB","MicrosporidiaDB", "PiroplasmaDB")
plasmodb = c("PadleriG01","PbergheiANKA","PbillcollinsiG01","PblacklockiG01","Pchabaudichabaudi","PcoatneyiHackeri","PcynomolgiB","PcynomolgiM","Pfalciparum3D7","Pfalciparum7G8","PfalciparumCD01","PfalciparumDd2","PfalciparumGA01","PfalciparumGB4","PfalciparumGN01","PfalciparumHB3","PfalciparumIT","PfalciparumKE01","PfalciparumKH01","PfalciparumKH02","PfalciparumML01","PfalciparumSD01","PfalciparumSN01","PfalciparumTG01","PfragileNilgiri","PgaboniG01","PgaboniSY75","Pgallinaceum8A","PinuiSanAntonio1","PknowlesiH","PknowlesiMalayanPk1A","PmalariaeUG01","PovalecurtisiGH01","PpraefalciparumG01","PreichenowiCDC","PreichenowiG01","PrelictumSGS1-like","PvinckeipetteriCR","Pvinckeivinckeivinckei","PvivaxP01","PvivaxSal1","Pyoeliiyoelii17X","Pyoeliiyoelii17XNL","PyoeliiyoeliiYM")
cryptodb = c("Candersoni30847","Chominis30976","ChominisTU502","ChominisTU502_2012","ChominisUdeA01","CmeleagridisUKMEL1","CmurisRN66","CparvumIowaII","CtyzzeriUGA55", "Cubiquitum39726","CveliaCCMP2878", "GniphandrodesUnknown", "VbrassicaformisCCMP3155")
giardiadb = c("GintestinalisAssemblageADH", "GintestinalisAssemblageAWB", "GintestinalisAssemblageBGS", "GintestinalisAssemblageBGS_B", "GintestinalisAssemblageEP15", "SsalmonicidaATCC50377") # No AssemblageA
tritrypdb = c("BayalaiB08-376","CfasciculataCfCl","EmonterogeiiLV88","LaethiopicaL147", "LarabicaLEM1108", "LbraziliensisMHOMBR75M2903", "LbraziliensisMHOMBR75M2904", "LdonovaniBPK282A1", "LenriettiiLEM3045", "LgerbilliLEM452","LinfantumJPCM5", "LmajorFriedlin", "LmajorLV39c5", "LmajorSD75.1", "LmexicanaMHOMGT2001U1103", "LpanamensisMHOMCOL81L13","LpanamensisMHOMPA94PSC1", "LpyrrhocorisH10", "LseymouriATCC30220", "LspMARLEM2494", "LtarentolaeParrotTarII", "LtropicaL590", "LturanicaLEM423", "PconfusumCUL13","TbruceigambienseDAL972", "TbruceiLister427", "TbruceiTREU927", "TcongolenseIL3000", "TcruziCLBrener", "TcruziCLBrenerEsmeraldo-like", "TcruziCLBrenerNon-Esmeraldo-like", "TcruzicruziDm28c","TcruziDm28c", "TcruzimarinkelleiB7", "TcruziSylvioX10-1", "TcruziSylvioX10-1-2012","TevansiSTIB805", "TgrayiANR4", "TrangeliSC58", "TvivaxY486", "TtheileriEdinburgh")
# MUST ASK PERMISSION TO USE MANY OF THE TRITRYP GENOMES - S.M. Beverley at Wash U
trichdb = c("TvaginalisG3")
amoebadb =  c("AcastellaniiNeff", "EdisparSAW760", "EhistolyticaHM1IMSS-A", "EhistolyticaHM1IMSS-B", "EhistolyticaHM1IMSS", "EhistolyticaHM3IMSS", "EhistolyticaKU27", "EinvadensIP1", "EmoshkovskiiLaredo", "EnuttalliP19", "NfowleriATCC30863") # ALL
toxodb = c("CcayetanensisCHN_HEN01", "CsuisWienI","EacervulinaHoughton", "EbrunettiHoughton", "EfalciformisBayerHaberkorn1970", "EmaximaWeybridge", "EmitisHoughton", "EnecatrixHoughton", "EpraecoxHoughton", "EtenellaHoughton", "HhammondiHH34", "NcaninumLIV", "SneuronaSN3", "SneuronaSOSN1", "TgondiiARI", "TgondiiFOU", "TgondiiGAB2-2007-GAL-DOM2", "TgondiiGT1", "TgondiiMAS", "TgondiiME49", "Tgondiip89", "TgondiiRH", "TgondiiRUB", "TgondiiTgCatPRC2", "TgondiiVAND", "TgondiiVEG")
microsporidiadb = c("AalgeraePRA109", "AalgeraePRA339", "AspWSBS2006","EaedisUSNM41457", "EbieneusiH348", "EcanceriGB1","EcuniculiEC1", "EcuniculiEC2", "EcuniculiEC3","EcuniculiEcunIII-L","EcuniculiGBM1", "EhellemATCC50504", "EhellemSwiss", "EhepatopenaeiTH1","EintestinalisATCC50506", "EromaleaeSJ2008","Heriocheircanceri","HeriocheirGB1", "MdaphniaeUGP3", "NausubeliERTm2", "NausubeliERTm6", "NbombycisCQ1", "NceranaeBRL01", "NceranaePA08_1199","NdisplodereJUm2807","NparisiiERTm1", "NparisiiERTm3", "OcolligataOC4", "PneurophiliaMK1", "Slophii42_110", "ThominisUnknown", "VcorneaeATCC50505", "Vculicisfloridensis")
piroplasmadb = c("BbigeminaBOND", "BbovisT2Bo", "BmicrotiRI","BovataMiyake", "CfelisWinnie", "TannulataAnkara", "TequiWA", "TorientalisShintoku", "TparvaMuguga")
fungidb = c("AaculeatusATCC16872", "AbrasiliensisCBS101740", "AcampestrisIBT28561", "Acandida2VRR", "AcarbonariusITEM5010", "AclavatusNRRL1", "AfischeriNRRL181","AflavusNRRL3357","AfumigatusA1163","AfumigatusAf293", "AglaucusCBS516.65","AinvadansNJM9701", "AlaibachiiNc14", "AluchuensisCBS106.47", "AmacrogynusATCC38327", "AnidulansFGSCA4", "AnigerATCC1015", "AnigerCBS513-88", "AnovofumigatusIBT16806", "AochraceoroseusIBT24754","AoryzaeRIB40", "AsteyniiIBT23096", "AsydowiiCBS593.65", "AterreusNIH2624", "AtubingensisCBS134.48", "AversicolorCBS583.65", "AwentiiDTO134E9","AzonataCBS506.65", "BcinereaB05-10", "BdendrobatidisJEL423", "CalbicansSC5314", "CalbicansSC5314_B", "CalbicansWO1", "CaurisB8441", "Ccinereaokay7-130", "CdeuterogattiiR265", "CgattiiCA1873", "CgattiiEJB2", "CgattiiIND107", "CgattiiWM276", "CglabrataCBS138", "CimmitisH538-4", "CimmitisRS", "ClusitaniaeATCC42720", "CneoformansB-3501A", "CneoformansH99", "CneoformansJEC21", "CneoformansKN99", "CposadasiiC735deltSOWgp", "CposadasiiRMSCC3488", "CposadasiiRMSCC3700", "CposadasiiSilveira", "FfujikuroiIMI58289", "FgraminearumPH-1", "Foxysporum26406", "Foxysporum4287", "Foxysporum54006", "FoxysporumFo47", "Foxysporumrace1", "Foxysporumrace4", "Fverticillioides7600", "HarabidopsidisEmoy2", "HcapsulatumG186AR", "HcapsulatumG217B", "HcapsulatumH143", "HcapsulatumH88", "HcapsulatumNAm1", "McircinelloidesCBS277-49", "MglobosaCBS7966", "Mlarici-populina98AG31", "Moryzae70-15", "MoryzaeBR32", "NcrassaOR74A", "NdiscretaFGSC8579", "NtetraspermaFGSC2508", "PaphanidermatumDAOMBR444", "ParrhenomanesATCC12531", "PblakesleeanusNRRL1555", "PbrasiliensisPb03", "PbrasiliensisPb18", "PcapsiciLT1534", "PchrysosporiumRP-78", "PcinnamomiCBS144-22", "PgraminisCRL75-36-700-3", "PinfestansT30-4", "PirregulareDAOMBR486", "PiwayamaiDAOMBR242034", "PjiroveciiSE8", "PlutziiPb01", "PparasiticaINRA-310", "PramorumPr-102", "PrubensWisconsin54-1255", "PsojaeP6497", "PultimumBR650", "PultimumDAOMBR144", "PvexansDAOMBR484", "RdelemarRA99-880", "ScerevisiaeS288c", "SdiclinaVS20", "SjaponicusyFS275", "Smacrosporak-hell", "SoctosporusyFS286", "SparasiticaCBS223", "Spombe972h", "SpunctatusDAOMBR117", "SreilianumSRZ2", "Sschenckii1099-18", "Ssclerotiorum1980UF-70", "TmarneffeiATCC18224", "TmesentericaDSM1558", "TreeseiQM6a", "TstipitatusATCC10500", "Umaydis521", "Ureesii1704", "YlipolyticaCLIB122", "ZtriticiIPO323",
    "AkawachiiIFO4308", "AnigerN402ATCC64974","CgattiiNT10","CparapsilosisCDC317","FproliferatumET1", "LprolificansJHH5317", "SapiospermumIHEM14462","Sbrasiliensis5110","YlipolyticaCLIB89W29")

# annotated proteins
genomes = vector(mode = "list", length = 9)
genomes[[1]] = plasmodb; genomes[[2]] = cryptodb; genomes[[3]] = giardiadb; genomes[[4]] =  tritrypdb; genomes[[5]] = trichdb; genomes[[6]] = amoebadb; genomes[[7]] = toxodb; genomes[[8]] = microsporidiadb; genomes[[9]] = piroplasmadb# ; genomes[[10]] = fungidb
# setwd("/scratch/glm5uh/eupathdb2models/genomes_from_EuPathDB/version40/annotated")
# wd_String = "/Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/data"
wd_String = "/home/mac9jc/paradigm/data"
# stop_counter = 0
for (i in 1:length(genomes)) { genomes_list = genomes[[i]]; db_name_uc = DataBase_uc[i]; db_name_lc = database_lc[i]
  for (x in genomes_list) {
    release_string = paste(paste(".org/common/downloads/release-",version, sep = ""),"/", sep = "")
    upstream_link = paste(paste("http://",db_name_lc,sep = ""),release_string, sep = "")
    middle_link = paste(paste(x,'/fasta/data/', sep = ""),db_name_uc, sep = "") 
    downstream_link = paste(paste(version, x, sep = "_"),"_AnnotatedProteins.fasta", sep = "")
    link = (paste(paste(upstream_link,middle_link,sep = ""), downstream_link,sep = "-")) 
    print(link)
    destination = paste(wd_String,paste(x,paste("annotated_Feb2018",".fasta", sep = ""), sep = "_"), sep = "/")
    print(destination)
    if (!url.exists(link)) {
        print('------------------------------------------------------')
        print('DOES NOT EXIST')
        print(link)
        list_o_broken_links = c(list_o_broken_links, link)}
    else {download.file(url = link, destfile = destination, method = "wget", quiet = T)}
  #  stop_counter = 1 + stop_counter
 #   if (stop_counter > 5) {break}
  }
#  if (stop_counter > 5) {break}
}
