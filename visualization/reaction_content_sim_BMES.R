
# # calc # of reaction shared by all models
# # cals most similar and most different models
# #### pca of reaction presence matrix, color code by genus

### USE GAPFILLED VERSION for paper
df = read.csv("/Users/maureencarey/local_documents/work/comparative_parasite_models/paradigm/data/reaction_matrix_july25.csv")
rownames(df) = df$X; df$X = NULL
df = df[,-which(colnames(df) == 'TgondiiRH')]
num_rxns_in_all = dim(subset(df, rowSums(df)==162))[1]
num_rxns_in_most = dim(subset(df, rowSums(df)>80))[1]
percent_in_most = num_rxns_in_most/dim(df)[1]

# d = dist(df)
# hc = hclust(d)
# fit <- cmdscale(d,eig=TRUE) # k is the number of dim
# x <- fit$points[,1]; y <- fit$points[,2]; colnames(fit$points) = c('x', 'y')
# rxn_occurance = ggplot(data = as.data.frame(fit$points), aes(x = x, y = y)) + geom_point() + xlab('Coordinate 1') + ylab('Coordinate 2')

d = dist(t(df))
df = as.data.frame(as.matrix(d))
df$species1 <- rownames(df)
melted_df = melt(df, id.vars = c('species1'));
melted_df2 = subset(melted_df, value<5);
melted_df1 = subset(melted_df2, value>0);

melted_df0 = subset(melted_df1, value<2);
melted_max = subset(melted_df, value>50)
# hc = hclust(d)
fit <- cmdscale(d,eig=TRUE) # k is the number of dim
x <- fit$points[,1]; y <- fit$points[,2]; colnames(fit$points) = c('x', 'y')


approx_genus= substr(rownames(fit$points),1,1)

rbc_para = ifelse(approx_genus == 'P',1,0) #"PneurophiliaMK1""PconfusumCUL13"  included
rbc_para[rownames(fit$points)=='PneurophiliaMK1'] = 0
rbc_para[rownames(fit$points)=='PconfusumCUL13'] = 0

crypto = ifelse(approx_genus == 'C',1,0)
crypto[rownames(fit$points)=='CfelisWinnie'] = 0
crypto[rownames(fit$points)=='CsuisWienI'] = 0
crypto[rownames(fit$points)=='CveliaCCMP2878'] = 0
crypto[rownames(fit$points)=='CcayetanensisCHN_HEN01'] = 0
crypto[rownames(fit$points)=='CfasciculataCfCl'] = 0

toxo = ifelse(substring(rownames(fit$points),1,7) == 'Tgondii' ,1, 0)

genus_color = ifelse(rbc_para == 1,"red",
                     ifelse(toxo == 1, "burlywood4",
                            ifelse(crypto == 1,"orange","gray")))

approx_genus2 = ifelse(toxo == 1, "toxo", approx_genus)
extracellular_para = ifelse(approx_genus == 'L',1,ifelse(approx_genus2 == 'T',1,0))
#epithelial_para

p0 = ggplot(data = as.data.frame(fit$points), aes(x = x, y = y)) + 
  geom_point(aes(color = as.factor(genus_color)), alpha = .5) + 
  scale_colour_manual(values = levels(as.factor(genus_color))) +
  xlab('Coordinate 1') + ylab('Coordinate 2') + guides(color = FALSE) + 
  ggtitle('Color by genus')
p1 = ggplot(data = as.data.frame(fit$points), aes(x = x, y = y)) + 
  geom_point(aes(color = as.factor(extracellular_para)), alpha = 0.5) + 
  scale_colour_manual(values = c("gray","blue")) +
  xlab('Coordinate 1') + ylab('Coordinate 2') + guides(color = FALSE, alpha = FALSE) + 
  ggtitle('Extracellular parasites in blue')
# p2 = ggplot(data = as.data.frame(fit$points), aes(x = x, y = y)) + geom_point(aes(color = as.factor(epithelial_para))) + xlab('Coordinate 1') + ylab('Coordinate 2') + guides(color = FALSE) + ggtitle('Color by epithelial cell invading parasite')
p3 = ggplot(data = as.data.frame(fit$points), aes(x = x, y = y)) + 
  geom_point(aes(color = as.factor(rbc_para)), alpha = 0.5) + 
  scale_colour_manual(values = c("gray","red")) +
  xlab('Coordinate 1') + ylab('Coordinate 2') + guides(color = FALSE, alpha = FALSE)  + 
  ggtitle('RBC-invading parasites in red')

p_final = ggarrange(p0, p1, p3,
                    labels = c("A", "B", "C"),
                    ncol = 3, nrow = 1)
ggsave("~/local_documents/work/dissertation/dissertation/figures/reaction_content_similarity_alpha.pdf", p_final, width = 10.5, height = 3, units = "in", dpi = 600)
#

crypto = ifelse(approx_genus2 == 'E',0,1)

ggplot(data = as.data.frame(fit$points), aes(x = x, y = y)) + 
  geom_point(aes(color = as.factor(crypto)), alpha = .5) + 
  xlab('Coordinate 1') + ylab('Coordinate 2') + guides(color = FALSE) + 
  ggtitle('Color by genus')
