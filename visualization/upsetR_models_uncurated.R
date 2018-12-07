library(UpSetR)
library(tidyverse)
library(seqinr)
library(Biostrings)
#df = read.csv("/Users/maureencarey/local_documents/work/comparative_parasite_models/data/models/reaction_matrix.csv")
#df = read.csv("/Users/maureencarey/local_documents/work/comparative_parasite_models/code/data/genomes_from_EuPathDB/version37/reaction_matrix_pruned.csv")
df = read.csv("/Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/data/reaction_matrix_july25.csv")

# 
# dist_df = df
# rownames(dist_df) = df$X; dist_df$X = NULL
# colnames(dist_df) = sub("_May.*", "", colnames(dist_df))
# dist_df = dist(t(dist_df))
# hc = hclust(dist_df)
# plot(hc, hang = -1, cex = 0.5)
# library("ape")
# plot(as.phylo(hc),  cex = 0.4)
# plot(as.phylo(hc), type = "fan")

path = "~/local_documents/work/comparative_parasite_models/paradigm/data/genomes_from_EuPathDB_version38/"
file.names <- dir(path, pattern ="2018.fasta")
num_genes = data.frame(species=rep(NA,length(file.names)), number_AA=rep(NA,length(file.names)))
for(i in 1:length(file.names)){
  num_genes$species[i] = file.names[i]
  fastaFile <- readAAStringSet(paste(path,file.names[i],sep = ""))
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  num_genes$number_AA[i] = length(seq_name)}

#colnames(df) = gsub("_May2018.xml","",colnames(df))
#colnames(df) = gsub("_May20.tsv","",colnames(df))
num_genes$species = gsub("_Jul2018.fasta","",num_genes$species)
rownames(df) = df$X; df$X = NULL

plasmodb = c("PadleriG01","PbergheiANKA","PbillcollinsiG01","PblacklockiG01","Pchabaudichabaudi","PcoatneyiHackeri","PcynomolgiB","PcynomolgiM","Pfalciparum3D7","PfalciparumIT","PfragileNilgiri","PgaboniG01","PgaboniSY75","Pgallinaceum8A","PinuiSanAntonio1","PknowlesiH","PknowlesiMalayanPk1A","PmalariaeUG01","PovalecurtisiGH01","PpraefalciparumG01","PreichenowiCDC","PreichenowiG01","PrelictumSGS1.like","PvinckeipetteriCR","Pvinckeivinckeivinckei","PvivaxP01","PvivaxSal1","Pyoeliiyoelii17X","Pyoeliiyoelii17XNL","PyoeliiyoeliiYM") # ALL
cryptodb = c("Candersoni30847","Chominis30976","ChominisTU502","ChominisTU502_2012","ChominisUdeA01","CmeleagridisUKMEL1","CmurisRN66","CparvumIowaII","CtyzzeriUGA55","Cubiquitum39726","CveliaCCMP2878", "GniphandrodesUnknown", "VbrassicaformisCCMP3155")
# removed Cbaileyi and Chominis37999 ChominisUKH1 since they don't have an annotated genome yet
giardiadb = c("GintestinalisAssemblageADH", "GintestinalisAssemblageAWB", "GintestinalisAssemblageBGS", "GintestinalisAssemblageBGS_B", "GintestinalisAssemblageEP15", "SsalmonicidaATCC50377") # No AssemblageA
tritrypdb = c("BayalaiB08.376","CfasciculataCfCl","EmonterogeiiLV88","LaethiopicaL147", "LarabicaLEM1108", 
              "LbraziliensisMHOMBR75M2903", "LbraziliensisMHOMBR75M2904", "LdonovaniBPK282A1", 
              "LenriettiiLEM3045", "LgerbilliLEM452","LinfantumJPCM5", "LmajorFriedlin", "LmajorLV39c5", 
              "LmajorSD75", "LmexicanaMHOMGT2001U1103", "LpanamensisMHOMCOL81L13","LpanamensisMHOMPA94PSC1", 
              "LpyrrhocorisH10", "LseymouriATCC30220", "LspMARLEM2494", "LtarentolaeParrotTarII", 
              "LtropicaL590", "LturanicaLEM423", "PconfusumCUL13","TbruceigambienseDAL972", "TbruceiLister427",
              "TbruceiTREU927", "TcongolenseIL3000", "TcruziCLBrener", "TcruziCLBrenerEsmeraldo.like", 
              "TcruziCLBrenerNon.Esmeraldo.like", "TcruzicruziDm28c","TcruziDm28c", "TcruzimarinkelleiB7", 
              "TcruziSylvioX10.1", "TcruziSylvioX10.1.2012","TevansiSTIB805", "TgrayiANR4", "TrangeliSC58", 
              "TvivaxY486", "TtheileriEdinburgh")
# remove LamazonensisMHOMBR71973M2269 LdonovaniBHU1220 TcruziEsmeraldo TcruziJRcl4 TcruziTulacl2
trichdb = c("TvaginalisG3")
amoebadb =  c("AcastellaniiNeff", "EdisparSAW760", "EhistolyticaHM1IMSS.A", "EhistolyticaHM1IMSS.B", "EhistolyticaHM1IMSS", "EhistolyticaHM3IMSS", "EhistolyticaKU27", "EinvadensIP1", "EmoshkovskiiLaredo", "EnuttalliP19", "NfowleriATCC30863") # ALL
# removed AcastellaniiMa AculbertsoniA1 AlenticulataPD2S AlugdunensisL3a Amauritaniensis1652 ApalestinensisReich AquinaVil3 ArhysodesSingh AspGalka AspIncertaesedis AspT4b-type AtriangularisSH621 EhistolyticaDS4 EhistolyticaHM1CA EhistolyticaKU48 EhistolyticaKU50 EhistolyticaMS96 AtriangularisSH621 EhistolyticaDS4 AquinaVil3    AastronyxisUnknown EhistolyticaRahman
toxodb = c("CcayetanensisCHN_HEN01", "CsuisWienI","EacervulinaHoughton", "EbrunettiHoughton", "EfalciformisBayerHaberkorn1970", "EmaximaWeybridge", "EmitisHoughton", "EnecatrixHoughton", "EpraecoxHoughton", "EtenellaHoughton", "HhammondiHH34", "NcaninumLIV", "SneuronaSN3", "SneuronaSOSN1", "TgondiiARI", "TgondiiFOU", "TgondiiGAB2.2007.GAL.DOM2", "TgondiiGT1", "TgondiiMAS", "TgondiiME49", "Tgondiip89", "TgondiiRH", "TgondiiRUB", "TgondiiTgCatPRC2", "TgondiiVAND", "TgondiiVEG")
# removed TgondiiCAST TgondiiCOUG TgondiiCtCo5 TgondiiTgCATBr5 TgondiiTgCATBr9 TgondiiTgCkUg2   
microsporidiadb = c("AalgeraePRA109", "AalgeraePRA339", "EaedisUSNM41457", "EbieneusiH348", "EcuniculiEC1", "EcuniculiEC2", "EcuniculiEC3", "EcuniculiGBM1", "EhellemATCC50504", "EhellemSwiss", "EintestinalisATCC50506","EromaleaeSJ2008","MdaphniaeUGP3", "NausubeliERTm2", "NausubeliERTm6", "NbombycisCQ1", "NceranaeBRL01", "NdisplodereJUm2807","NparisiiERTm1", "NparisiiERTm3", "OcolligataOC4", "PneurophiliaMK1", "Slophii42_110", "ThominisUnknown", "VcorneaeATCC50505", "Vculicisfloridensis")
# removed AalgeraeUndeen HtvaerminnensisOER-3-3
piroplasmadb = c("BbigeminaBOND", "BbovisT2Bo", "BmicrotiRI","BovataMiyake", "CfelisWinnie", "TannulataAnkara", "TequiWA", "TorientalisShintoku", "TparvaMuguga") # ALL



database = NA
for (i in 1:ncol(df)) {
  
  if (colnames(df)[i]%in%plasmodb) { y = 'PlasmoDB'
  } else if (colnames(df)[i]%in%cryptodb) { y = 'CryptoDB'
  } else if (colnames(df)[i]%in%giardiadb) { y = 'GiardiaDB'
  } else if (colnames(df)[i]%in%tritrypdb) { y = 'TriTrypDB'
  } else if (colnames(df)[i]%in%trichdb) { y = 'TrichDB'
  } else if (colnames(df)[i]%in%amoebadb) { y = 'AmoebaDB'
  } else if (colnames(df)[i]%in%toxodb) { y = 'ToxoDB'
  } else if (colnames(df)[i]%in%microsporidiadb) { y = 'MicrosporidiaDB'
  } else if (colnames(df)[i]%in%piroplasmadb) { y = 'PiroplasmaDB'
  } else { y = 'error'
  print(colnames(df)[i])}
  database[i] = y
}

library(reshape2)
cn = colnames(df)
metadata = as.data.frame(cbind(cn,database))
colnames(metadata) = c('sets','database')
metadata[] <- lapply(metadata, as.character)
for (i in 1:length(metadata$sets)) {
  n <- 2
  x = metadata$sets[i]
  metadata$sets[i] = gsub("([[:lower:]])([[:upper:]])", "\\1 \\2", paste(substr(x, 1, n-1), ". ", substr(x, n, nchar(x)), sep = ""))
}
for (i in 1:length(colnames(df))) {
  n <- 2
  x = colnames(df)[i]
  colnames(df)[i] = gsub("([[:lower:]])([[:upper:]])", "\\1 \\2", paste(substr(x, 1, n-1), ". ", substr(x, n, nchar(x)), sep = ""))
}
for (i in 1:length(num_genes$species)) {
  n <- 2
  x = num_genes$species[i]
  num_genes$species[i] = gsub("([[:lower:]])([[:upper:]])", "\\1 \\2", paste(substr(x, 1, n-1), ". ", substr(x, n, nchar(x)), sep = ""))
}

colnames(df)[colnames(df)=='L. major SD75'] = 'L. major SD75.1'
metadata$sets[metadata$sets == 'L. major SD75'] = 'L. major SD75.1'
num_genes$species = gsub("-", ".", num_genes$species)
#num_genes = num_genes[!(num_genes$species == "T. cruzi CLBrener" | num_genes$species == "T. gondii RH"),]
if (length(setdiff(num_genes$species,colnames(df)))> 0) {print('difference')}
if (length(setdiff(num_genes$species,metadata$sets))> 0) {print('difference')}

md = metadata[order(metadata$database, metadata$sets),]
df2 = df[,md$sets]


contains_reactions = data.frame(species = colnames(df2), contains_rxns = rep(NA, length(colnames(df2))))
for (i in 1:length(colnames(df2))) {
  contains_reactions$contains_rxns[i] = list(rownames(df2)[df2[,i] != 0])
}
# setdiff(c(1,2,3),c(2,3,4)) = 1
unique_reactions = data.frame(species = colnames(df2), unique_rxns = rep(NA, length(colnames(df2))))
for (i in 1:length(colnames(df2))) {
  unique_reactions$unique_rxns[i] = list(setdiff(unlist(contains_reactions$contains_rxns[i]), unlist(contains_reactions$contains_rxns[-i])))
}

num_reactions = data.frame(species = colnames(df2), num_reactions = colSums(df2))

colnames(md) = c('species','database')
df = merge(num_reactions,num_genes, by = 'species')
df = merge(df, unique_reactions, by = 'species')
df = merge(df, md, by = 'species')
num_unique_rxns = rep(NA, nrow(df))
df = cbind(df, num_unique_rxns)
for (i in 1:nrow(df)) {
  if (length(df$unique_rxns[i]) == 0) { df$num_unique_rxns[i] = 0 
  } else {df$num_unique_rxns[i]  = str_count(df$unique_rxns[i], ',') + 1} }

df$species <- factor(df$species, levels = df$species[order(df$database, df$species)])

p = ggplot() +
  geom_point(data = df, aes(x = species, y = number_AA, fill = as.factor(database)), size = 2) +  
  guides(fill = F) +
  geom_bar(data = df, aes(x = species, y = num_reactions, fill = as.factor(database)), color = "black", size = .1, stat = "identity") + 
  geom_point(data = df, aes(x = species, y = num_unique_rxns), color = "black", size = 1) + 
  scale_y_log10(minor_breaks = seq(0, 10000, 1000),
                breaks = c(1.0,10,100,1000,10000),
                labels = function(x) format(x, scientific = FALSE), position = "right") + coord_flip() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.background = element_blank(),
        panel.border= element_rect(color = "black", fill = NA, size=1),
        panel.grid.major.x = element_line(color = "grey"),
        panel.grid.minor.x = element_line(color = "grey"),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.spacing.x = unit(2, "lines"),
        legend.title=element_blank(),
        legend.position = c(.82, 0.75),
        legend.text=element_text(size=8)) 
#setwd(paste(base_directory,'/figures',sep=""))
#ggsave("temp3.pdf", p, width = 9, height = 30, units = "in", dpi = 600)


rects = df[order(df$database),c(1,5)]
c1 = sort((df$database))
xstart = 1:length(c1) - 0.5
xend = 1:length(c1) + 0.5
c_df = data.frame(as.character(c(t(c1))), as.numeric(xstart), as.numeric(xend))
colnames(c_df) = c('database','xstart','xend')

p = ggplot()  + 
  geom_point(data = df, aes(x = species, y = number_AA, fill = as.factor(database)), size = 1, shape = 2) +
  geom_rect(data = c_df, aes(ymin=1, ymax=350000, xmin=xstart,
                             xmax=xend, fill=database), alpha =0.5) +
  geom_point(data = df, aes(x = species, y = number_AA, fill = as.factor(database)), size = 1, shape = 2) +
  guides(fill = F) +
  geom_bar(data = df, aes(x = species, y = num_unique_rxns, fill = as.factor(database)), color = "black", size = .1, stat = "identity") + 
  geom_point(data = df, aes(x = species, y = num_reactions), color = "black", size = 1, shape = 1) + 
  scale_y_log10(minor_breaks = c(5, 50,500,5000,50000),
                breaks = c(0,1,10,100,1000,10000),
                labels = function(x) format(x, scientific = FALSE), position = "right") + coord_flip() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        panel.background = element_blank(),
        panel.border= element_rect(color = "black", fill = NA, size=1),
        panel.grid.major.x = element_line(color = "grey"),
        panel.grid.minor.x = element_line(color = "grey"),
        panel.grid.major.y = element_blank(), 
        panel.grid.minor.y = element_blank(),
        panel.spacing.x = unit(2, "lines"),
        legend.title=element_blank(),
        legend.position = c(.82, 0.75),
        legend.text=element_text(size=8)) 

#setwd(paste(base_directory,'/figures',sep=""))
# ggsave("temp4.pdf", p, width = 9, height = 30, units = "in", dpi = 600)  

p2 = p + theme(axis.title.x = element_blank(),
               axis.title.y = element_blank(),
               axis.text.x = element_text(size = 6),
               axis.text.y = element_text(size = 6),
               panel.background = element_blank(),
               panel.border= element_rect(color = "black", fill = NA, size=0.5),
               panel.grid.major.x = element_line(color = "grey"),
               panel.grid.minor.x = element_line(color = "grey"),
               panel.grid.major.y = element_blank(), 
               panel.grid.minor.y = element_blank(),
               panel.spacing.x = unit(2, "lines"),
               legend.title=element_blank(),
               legend.position = c(.82, 0.75),
               legend.text=element_text(size=7)) 
p2poster = ggplot() + 
  geom_point(data = df, aes(x = species, y = number_AA, fill = as.factor(database)), size = 3, shape = 2) +
  geom_rect(data = c_df, aes(ymin=1, ymax=350000, xmin=xstart, xmax=xend, fill=database), alpha =0.5) +
  geom_point(data = df, aes(x = species, y = number_AA, fill = as.factor(database)), size = 3, shape = 2) +
  guides(fill = F) +
  geom_bar(data = df, aes(x = species, y = num_unique_rxns, fill = as.factor(database)), color = "black", size = .1, stat = "identity") + 
  geom_point(data = df, aes(x = species, y = num_reactions), color = "black", size = 3, shape = 1) + 
  geom_point(data = df, aes(x = species, y = num_reactions), color = "black", size = 3, shape = 1) + 
  scale_y_log10(minor_breaks = c(5, 50,500,5000,50000),
                breaks = c(0,1,10,100,1000,10000),
                labels = function(x) format(x, scientific = FALSE), position = "right") + coord_flip() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_blank(),
        panel.background = element_blank(),
        panel.border= element_rect(color = "black", fill = NA, size=0.5),
        panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(),
        panel.spacing.x = unit(2, "lines"), legend.title=element_blank()) 
ggsave("~/local_documents/work/dissertation/dissertation/figures/model_overview_poster.pdf", p2poster, width = 8, height = 29.9, units = "in", dpi = 600)              

ggsave("~/local_documents/work/dissertation/dissertation/figures/model_overview_by_database25.pdf", p2, width = 3, height = 11.5, units = "in", dpi = 600)              

df = df[which(df$species != 'T. gondii RH'),]
mod <- lm(num_reactions ~ log10(number_AA), data = df)
r = summary(mod)[['r.squared']]
pval = (summary(mod))$coefficients[2,4]
if (pval < 0.001) { pval = 'p < 0.001'}

lm_eqn = function(m) {
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  if (coef(m)[2] >= 0)  {
    eq <- substitute(~~italic(r)^2~"="~r2,l)
    #eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    #eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)     }
    eq <- substitute(~~italic(r)^2~"="~r2,l)     }
  as.character(as.expression(eq))}

g1 = ggplot(data = df, aes(x = number_AA, y = num_reactions)) + 
  geom_point() + 
  scale_x_log10(breaks = c(1.0,100, 1000,5000,10000,25000,50000), labels = function(x) format(x, scientific = FALSE)) +
  xlab('genome size\n(no. of ORFs)') + ylab('model size (no. of reactions)') +
  geom_smooth(method = 'lm', se = F) +
  annotate("text", x = 25000, y = 1070, size = 4, label = lm_eqn(mod), parse = T) +
  annotate("text", x = 25200, y = 850, size = 4, label = pval)

poster1 = ggplot(data = df, aes(x = number_AA, y = num_reactions)) + 
  geom_point() + 
  scale_x_log10(breaks = c(1.0,100, 1000,5000,10000,25000,50000), labels = function(x) format(x, scientific = FALSE)) +
  xlab('genome size\n(no. of ORFs)') + ylab('model size\n(no. of reactions)') +
  geom_smooth(method = 'lm', se = F) +
  annotate("text", x = 25000, y = 1170, size = 8, label = lm_eqn(mod), parse = T)+
  theme(axis.text = element_text(size = 18),
      axis.title=element_text(size=20))

g2 = ggplot() + geom_point(data = df, aes( x = number_AA, y = num_unique_rxns)) + 
  scale_x_log10(breaks = c(1.0,1000,10000,50000),
                labels = function(x) format(x, scientific = FALSE))+
  xlab('genome size\n(no. of ORFs)') + ylab('no. of unique reactions')

mod <- lm(num_reactions ~ num_unique_rxns, data = df)
r = summary(mod)[['r.squared']]
pval = (summary(mod))$coefficients[2,4]
if (pval < 0.001) { pval = 'p < 0.001'}

df[which(df$num_unique_rxns == max(df$num_unique_rxns)),]
df[which(df$number_AA == min(df$number_AA)),]

g3 = ggplot(data = df, aes( x = num_reactions, y = num_unique_rxns)) + geom_point()  + 
  xlab('model size\n(no. of reactions)') + ylab('no. of unique reactions') +
  geom_smooth(method = 'lm', se = F) +
  annotate("text", x = 1000, y = 55, size = 4, label = lm_eqn(mod), parse = T) +
  annotate("text", x = 1010, y = 45, size = 4, label = pval)
g3
poster = ggplot(data = df, aes( x = num_reactions, y = num_unique_rxns)) + geom_point()  + 
  xlab('model size\n(no. of reactions)') + ylab('no. of unique reactions')+
  theme(axis.text = element_text(size = 16),
        axis.title=element_text(size=16))

zeros = df[df$num_unique_rxns ==1,]
g4 = ggplot(data = zeros) + geom_histogram(aes(x = num_reactions), bins = 40, fill = "black")+ 
  xlab('model size\n(no. of reactions)') + ylab('no. of models with only\none unique reaction')

library(ggpubr)


# histogram of reaction freq (over some threshold maybe) and genome size
# reaction as point, frequency of in model v avg model size that it's in
# hypothesis less common reactions occur in bigger genomes?
df3 = read.csv("/Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/data/reaction_matrix_july25.csv")
rownames(df3) = df3$X
df3 = df3 %>% mutate(s = rowSums(.[-1]))
df3 = df3[,c(1,dim(df3)[2] )]

g5 = ggplot(data = df3, aes(x = s)) + 
  geom_histogram(bins = 160, fill = "black") +
  ylab('no. of reactions') + xlab('no. of models\ncontaining a reaction') + scale_y_sqrt()

poster5 = ggplot(data = df3, aes(x = s)) + 
  geom_histogram(bins = 160, fill = "black") +
  ylab('no. of reactions') + xlab('no. of models\ncontaining a reaction') + 
  theme(axis.text = element_text(size = 18),
        axis.title=element_text(size=20)) +
  scale_y_sqrt() 

p = ggarrange(g1, g3, g4, g5, # + rremove("x.text"), 
              labels = c("A", "B", "C", "D"),
              ncol = 2, nrow = 2)
ggsave("~/local_documents/work/dissertation/dissertation/figures/model_sum_comparison25.pdf", p, width = 8, height = 7, units = "in", dpi = 600)

p = ggarrange(poster1, poster5, # + rremove("x.text"), 
              labels = c("A", "B"),
              ncol = 1, nrow = 2)
ggsave("~/local_documents/work/dissertation/dissertation/figures/model_sum_poster.pdf", p, width = 6, height = 8, units = "in", dpi = 600)

