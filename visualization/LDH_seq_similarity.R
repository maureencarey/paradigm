library(msa)
library(Biostrings)
library(ggdendro)
library(seqinr)
library(ggplot2)
mySequences = readAAStringSet("~/local_documents/work/comparative_parasite_models/paradigm/data/LDH_eupathdb_july18.fasta.txt")
myAlignment = msa(mySequences)
myAlignment2 = msaConvert(myAlignment, type="seqinr::alignment")
d <- dist.alignment(myAlignment2, "identity")
d[is.na(d)] = 0
hc = hclust(d)

dend = as.dendrogram(hc); dend_data <- dendro_data(dend, type = "rectangle")

new_label = c()
for (i in 1:length(dend_data$labels$label)) {
  string = as.character(dend_data$labels$label[i])
  x = strsplit(string,' \\| ')[[1]]; new_label[i] = x[2] } #(paste(x[2],x[3],sep = ", "))}

col = unique(lapply(new_label, function(x) {strsplit(x,' ')[[1]][1]}))
col2 = palette(rainbow(length(col))); color_vec = c()
for (i in 1:length(new_label)) { genus = strsplit(new_label[i],' ')[[1]][1]
color_vec[i] = col2[which(col == genus)] }

p = ggplot(dend_data$segments) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend))+
  geom_text(data = dend_data$labels, aes(x, y, label = new_label, color = color_vec),
            hjust = 1, size = 1.3) + coord_flip() + ylim(-3,1)+ guides(color = F) +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank()) + xlab(NULL)+ ylab(NULL) +
  coord_flip()

ggsave("~/local_documents/work/dissertation/dissertation/figures/LDH_cluster.pdf", p, width = 3, height = 9, units = "in", dpi = 600)
