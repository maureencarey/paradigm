library(tidyverse)
library(ggpubr)
df = read_csv("/Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/data/ortho_annotations_per_genome_july23.csv"); df$X1 = NULL

df$species = substr(df$species,1,nchar(df$species)-21)
new_list = c()
for (i in 1:nrow(df)) { 
  x = str_sub(df$OrthoMCL_gene[i], 3, -3)
  new_list[[i]] = strsplit(x,"', '")}
newest_list = list()
for (i in 1:length(new_list)) {newest_list[[df$species[i]]] = new_list[i][[1]][[1]]}
rm(new_list,i,x)

plasmodb = c("PadleriG01","PbergheiANKA","PbillcollinsiG01","PblacklockiG01","Pchabaudichabaudi","PcoatneyiHackeri","PcynomolgiB","PcynomolgiM","Pfalciparum3D7","PfalciparumIT","PfragileNilgiri","PgaboniG01","PgaboniSY75","Pgallinaceum8A","PinuiSanAntonio1","PknowlesiH","PknowlesiMalayanPk1A","PmalariaeUG01","PovalecurtisiGH01","PpraefalciparumG01","PreichenowiCDC","PreichenowiG01","PrelictumSGS1-like","PvinckeipetteriCR","Pvinckeivinckeivinckei","PvivaxP01","PvivaxSal1","Pyoeliiyoelii17X","Pyoeliiyoelii17XNL","PyoeliiyoeliiYM") # ALL
cryptodb = c("Candersoni30847","Chominis30976","ChominisTU502","ChominisTU502_2012","ChominisUdeA01","CmeleagridisUKMEL1","CmurisRN66","CparvumIowaII","CtyzzeriUGA55","Cubiquitum39726","CveliaCCMP2878", "GniphandrodesUnknown", "VbrassicaformisCCMP3155")
# removed Cbaileyi and Chominis37999 ChominisUKH1 since they don't have an annotated genome yet
giardiadb = c("GintestinalisAssemblageADH", "GintestinalisAssemblageAWB", "GintestinalisAssemblageBGS", "GintestinalisAssemblageBGS_B", "GintestinalisAssemblageEP15", "SsalmonicidaATCC50377") # No AssemblageA
tritrypdb = c("BayalaiB08-376","CfasciculataCfCl","EmonterogeiiLV88","LaethiopicaL147", "LarabicaLEM1108", 
              "LbraziliensisMHOMBR75M2903", "LbraziliensisMHOMBR75M2904", "LdonovaniBPK282A1", "LenriettiiLEM3045", "LgerbilliLEM452","LinfantumJPCM5", "LmajorFriedlin", "LmajorLV39c5", 
              "LmajorSD75.1", "LmexicanaMHOMGT2001U1103", "LpanamensisMHOMCOL81L13","LpanamensisMHOMPA94PSC1", "LpyrrhocorisH10", "LseymouriATCC30220", "LspMARLEM2494", "LtarentolaeParrotTarII", 
              "LtropicaL590", "LturanicaLEM423", "PconfusumCUL13","TbruceigambienseDAL972", "TbruceiLister427", "TbruceiTREU927", "TcongolenseIL3000", "TcruziCLBrener", "TcruziCLBrenerEsmeraldo-like", 
              "TcruziCLBrenerNon-Esmeraldo-like", "TcruzicruziDm28c","TcruziDm28c", "TcruzimarinkelleiB7", "TcruziSylvioX10-1", "TcruziSylvioX10-1-2012","TevansiSTIB805", "TgrayiANR4", "TrangeliSC58", "TvivaxY486", "TtheileriEdinburgh")
# remove LamazonensisMHOMBR71973M2269 LdonovaniBHU1220 TcruziEsmeraldo TcruziJRcl4 TcruziTulacl2
trichdb = c("TvaginalisG3")
amoebadb =  c("AcastellaniiNeff", "EdisparSAW760", "EhistolyticaHM1IMSS-A", "EhistolyticaHM1IMSS-B", "EhistolyticaHM1IMSS", "EhistolyticaHM3IMSS", "EhistolyticaKU27", "EinvadensIP1", "EmoshkovskiiLaredo", "EnuttalliP19", "NfowleriATCC30863") # ALL
# removed AcastellaniiMa AculbertsoniA1 AlenticulataPD2S AlugdunensisL3a Amauritaniensis1652 ApalestinensisReich AquinaVil3 ArhysodesSingh AspGalka AspIncertaesedis AspT4b-type AtriangularisSH621 EhistolyticaDS4 EhistolyticaHM1CA EhistolyticaKU48 EhistolyticaKU50 EhistolyticaMS96 AtriangularisSH621 EhistolyticaDS4 AquinaVil3    AastronyxisUnknown EhistolyticaRahman
toxodb = c("CcayetanensisCHN_HEN01", "CsuisWienI","EacervulinaHoughton", "EbrunettiHoughton", "EfalciformisBayerHaberkorn1970", "EmaximaWeybridge", "EmitisHoughton", "EnecatrixHoughton", "EpraecoxHoughton", "EtenellaHoughton", "HhammondiHH34", "NcaninumLIV", "SneuronaSN3", "SneuronaSOSN1", "TgondiiARI", "TgondiiFOU", "TgondiiGAB2-2007-GAL-DOM2", "TgondiiGT1", "TgondiiMAS", "TgondiiME49", "Tgondiip89", "TgondiiRH", "TgondiiRUB", "TgondiiTgCatPRC2", "TgondiiVAND", "TgondiiVEG")
# removed TgondiiCAST TgondiiCOUG TgondiiCtCo5 TgondiiTgCATBr5 TgondiiTgCATBr9 TgondiiTgCkUg2   
microsporidiadb = c("AalgeraePRA109", "AalgeraePRA339", "EaedisUSNM41457", "EbieneusiH348", "EcuniculiEC1", "EcuniculiEC2", "EcuniculiEC3", "EcuniculiGBM1", "EhellemATCC50504", "EhellemSwiss", "EintestinalisATCC50506","EromaleaeSJ2008","MdaphniaeUGP3", "NausubeliERTm2", "NausubeliERTm6", "NbombycisCQ1", "NceranaeBRL01", "NdisplodereJUm2807","NparisiiERTm1", "NparisiiERTm3", "OcolligataOC4", "PneurophiliaMK1", "Slophii42_110", "ThominisUnknown", "VcorneaeATCC50505", "Vculicisfloridensis")
# removed AalgeraeUndeen HtvaerminnensisOER-3-3
piroplasmadb = c("BbigeminaBOND", "BbovisT2Bo", "BmicrotiRI","BovataMiyake", "CfelisWinnie", "TannulataAnkara", "TequiWA", "TorientalisShintoku", "TparvaMuguga") # ALL

sets = df$species
database = NA
for (i in 1:length(sets)) {
  if (sets[i]%in%plasmodb) { y = 'PlasmoDB'
  } else if (sets[i]%in%cryptodb) { y = 'CryptoDB'
  } else if (sets[i]%in%giardiadb) { y = 'GiardiaDB'
  } else if (sets[i]%in%tritrypdb) { y = 'TriTrypDB'
  } else if (sets[i]%in%trichdb) { y = 'TrichDB'
  } else if (sets[i]%in%amoebadb) { y = 'AmoebaDB'
  } else if (sets[i]%in%toxodb) { y = 'ToxoDB'
  } else if (sets[i]%in%microsporidiadb) { y = 'MicrosporidiaDB'
  } else if (sets[i]%in%piroplasmadb) { y = 'PiroplasmaDB'
  } else { y = 'error'
  print(sets[i])}
  database[i] = y }

rm(i,df)
library(UpSetR)
# upset(fromList(newest_list),nsets = 10, nintersects = 10, order.by = "freq", mainbar.y.label = "No. Genes Shared",
#       sets.x.label = "OrthoMCL annotations per genome", point.size = 0.3, line.size = 0.15, mb.ratio = c(0.2,0.8), show.numbers = "no", text.scale = 0.6)#, cutoff = 10000)

png("~/local_documents/work/dissertation/dissertation/figures/genome_sum2_A.png", width = 6, height = 2.1, units = "in", res = 1000)
upset(fromList(newest_list),nsets = 5, nintersects = 10, order.by = "freq", mainbar.y.label = "No. Gene\nAnnotations\nShared", mb.ratio = c(0.34,0.66),
                    sets.x.label = "OrthoMCL annotations per genome", point.size = 0.4, line.size = 0.25, show.numbers = "no", text.scale = 0.8)#, cutoff = 10000)
dev.off()

union = c()
#intersection = c()
for (i in 1:length(newest_list)) {
  genes_in_species = as.character(newest_list[[i]])
  species = names(newest_list)[[i]]
  print(species)
  if (species%in%plasmodb) { y = 'PlasmoDB'
  } else if (species%in%cryptodb) { y = 'CryptoDB'
  } else if (species%in%giardiadb) { y = 'GiardiaDB'
  } else if (species%in%tritrypdb) { y = 'TriTrypDB'
  } else if (species%in%trichdb) { y = 'TrichDB'
  } else if (species%in%amoebadb) { y = 'AmoebaDB'
  } else if (species%in%toxodb) { y = 'ToxoDB'
  } else if (species%in%microsporidiadb) { y = 'MicrosporidiaDB'
  } else if (species%in%piroplasmadb) { y = 'PiroplasmaDB'
  } else { y = 'error'}
  union[[y]] = c(union[[y]], genes_in_species) }
#   intersection[[y]] = c(intersection[[y]], genes_in_species) }

for (i in names(union)) { print(length(union[[i]]))}
for (i in names(union)) { union[[i]] = unique(union[[i]])}
for (i in names(union)) { print(length(union[[i]]))}

# for (i in names(intersection)) { print(length(intersection[[i]]))}
# for (i in 1:length(intersection)) { 
#   a = table(intersection[[i]])
#   m = max(a)
#   t = a[names(a)==m]
#   intersection[[i]] = 1 }
# for (i in names(intersection)) { print(length(intersection[[i]]))} # not working
# 
png("~/local_documents/work/dissertation/dissertation/figures/genome_sum2_B.png", width = 6, height = 3.5, units = "in", res = 1000)
upset(fromList(union),order.by = "freq", sets = unique(database), mainbar.y.label = "No. Gene\nAnnotations\nShared",
      sets.x.label = "OrthoMCL annotations per group", point.size = 0.4, line.size = 0.25, mb.ratio = c(0.2,0.8), show.numbers = "no", text.scale = 0.8) 
dev.off()
