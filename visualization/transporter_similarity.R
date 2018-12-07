library(ggplot2)
library(ggpubr)

# if made on desktop:
# transporter_dist = read.csv('/Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/data/transporter_similarity_before_gapfilling_july22.csv')
# rownames(transporter_dist) = transporter_dist$X; transporter_dist$X = NULL

# if made on rivanna
setwd("~/local_documents/work/comparative_parasite_models/paradigm/data/rivanna_results/")
# MUST ADD JUl2018_BiGG to LmajorSD75
temp = list.files(pattern="*Jul2018_BiGG.csv")
myfiles = lapply(temp, read.csv)
for (i in 1:length(myfiles)) {
  file = myfiles[[i]]; rownames(file) = file$X; file$X = NULL
  myfiles[[i]] = file}
library(tidyverse)
df = bind_rows(lapply(myfiles,rownames_to_column))
df = as.data.frame(df)
rownames(df) = df$rowname; df$rowname = NULL
df = df[ , order(names(df))]

df$DIY6_BayalaiB08.376_May20 = NULL
df$DIY6_ChominisTU502_May20 = NULL
df$DIY6_CveliaCCMP2878_May20 = NULL
df$DIY6_EbrunettiHoughton_May20 = NULL
df$DIY6_NausubeliERTm6_May20 = NULL
df$DIY6_PvivaxSal1_May20 = NULL
df$DIY6_TgondiiME49_May20 = NULL
df$DIY6_TgondiiRH_May20 = NULL
df$DIY6_TgondiiRH_Jul2018_BiGG = NULL
df = df[rownames(df) != 'DIY6_TgondiiRH_Jul2018_BiGG',]

for ( i in 1:nrow(df)) { for (j in 1:ncol(df)) {
    if (df[i,j] !=df[j,i] ){print(paste(paste('i=',i),paste('j=',j)))} } }

#plot(stats::hclust(as.dist(as.matrix(transporter_dist))))

# start here if on personal computer
for (j in 1:ncol(df)) {
  col_max = max(na.omit(df[,j]))
  for (i in 1:nrow(df)) {
    if (is.na(df[i,j])) { next }
    else {
      df[i,j] = df[i,j]/col_max }}}

fit <- cmdscale(as.dist(df),eig=TRUE) # k is the number of dim
x <- fit$points[,1]; y <- fit$points[,2]; colnames(fit$points) = c('x', 'y')
p = ggplot(data = as.data.frame(fit$points), aes(x = x, y = y)) + geom_point() + xlab('Coordinate 1') + ylab('Coordinate 2')

## FIX THIS APPROXIMATION
approx_genus = sapply(rownames(fit$points), function(x) {substring(x,6,6)})
extracellular_para = ifelse(approx_genus == 'L',1,ifelse(approx_genus == 'T',1,0))
extracellular_para = ifelse(substring(names(extracellular_para),3,9) == 'Tgondii' ,0, extracellular_para)
#epithelial_para = approx_genus
rbc_para = ifelse(approx_genus == 'P',1,0)

p0 = ggplot(data = as.data.frame(fit$points), aes(x = x, y = y)) + geom_point(aes(color = as.factor(approx_genus))) + xlab('Coordinate 1') + ylab('Coordinate 2') + guides(color = FALSE) + ggtitle('Color by genus')
p1 = ggplot(data = as.data.frame(fit$points), aes(x = x, y = y)) + geom_point(aes(color = as.factor(extracellular_para))) + xlab('Coordinate 1') + ylab('Coordinate 2') + guides(color = FALSE) + ggtitle('Color by intra/extracellular parasite')
# p2 = ggplot(data = as.data.frame(fit$points), aes(x = x, y = y)) + geom_point(aes(color = as.factor(epithelial_para))) + xlab('Coordinate 1') + ylab('Coordinate 2') + guides(color = FALSE) + ggtitle('Color by epithelial cell invading parasite')
p3 = ggplot(data = as.data.frame(fit$points), aes(x = x, y = y)) + geom_point(aes(color = as.factor(rbc_para))) + xlab('Coordinate 1') + ylab('Coordinate 2') + guides(color = FALSE)   + ggtitle('Color by RBC-invading parasite')

p_final = ggarrange(p0, p1, p3, 
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)
ggsave("~/local_documents/work/dissertation/dissertation/figures/transporter_similarity.pdf", p_final, width = 10.5, height = 3, units = "in", dpi = 600)
