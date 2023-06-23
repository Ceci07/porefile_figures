library(phyloseq)
library(ggplot2)
library(stringr)
library(dplyr)
library(Metrics)
library(viridis)

#esto esta aca /mnt/cive/csalazar/projects/3.porfile/real_data_public/matsuo/nanopore/filtered/samples
setwd("/mnt/raid2tb/bubble/tesis/capitulo2/Fig6/results/")

list.files()

c1 <- as.data.frame(read.table("COUNTS.tsv", sep = "\t", header = T))
t1 <- as.data.frame(read.table("TAXCLA.tsv", sep = "\t", header = T))
name <- as.data.frame(colnames(c1))
df <- as.data.frame(str_split_fixed(name$`colnames(c1)`, "_", 3))
df1 <- as.data.frame(cbind(name, df$V2))
colnames(df1) <- c("file", "id")
head(df1)

meta1 <- df1
meta1$round <- c('polish')
row.names(meta1) <- meta1$file
head(meta1)
dim(meta1)


c1 <- as.matrix(read.table("COUNTS.tsv", sep = "\t", header = T))
t1 <- as.matrix(read.table("TAXCLA.tsv", sep = "\t", header = T))

OTU1 = otu_table(c1, taxa_are_rows = TRUE)
TAX1 = tax_table(t1)
samples = sample_data(meta1)
OTU1
TAX1

ps1 = phyloseq(OTU1, TAX1, samples)
ps1

ps <- psmelt(ps1)

pp <- ps[which(ps$id == "DRR225046"),]

top20 <- names(sort(taxa_sums(ps1), decreasing=TRUE))[1:10000]
ps.top20 <- transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

graph1 <- psmelt(ps.top20)
head(graph1)
tab2 <- graph1
tab2$step <- c("match")

write.table(tab2, "matched_abundance_tab.tsv", sep = "\t", row.names = F, quote = F)

mock <- graph1[which(graph1$id == "DRR225046"),]
unique(mock$Species)


p <- ggplot(mock, aes(x = id, y = Abundance, fill = Species))
m1 <- p +  geom_bar(aes(), stat="identity", position="fill") + 
  # scale_fill_manual(values = coul) +
  theme(legend.position="right", axis.line = element_line(colour = "black"), 
        text = element_text(size = 10), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 10, angle = 90),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  guides(fill=guide_legend(ncol = 5)) + ggtitle("porefile") +
  xlab("") + ylab("Relative Abundance") + labs(fill = "Species") 
m1


target <- c("Cereibacter sphaeroides", "Streptococcus mutans", "Enterococcus faecalis", "Bacillus cereus", "Clostridium beijerinckii",
            "Bifidobacterium adolescentis", "Escherichia coli", "Lactobacillus gasseri", "Staphylococcus epidermidis", "Streptococcus mutans", "Deinococcus radiodurans")

mock2 <- filter(mock, Species %in% target)
mock2 <- mock2 %>% group_by(Species) %>%
  summarise(sum = sum(Abundance))
max(mock2$sum)
min(mock2$sum)
unique(mock2$Species)

mock2$expected <- c("0.10")
mock2 <- as.data.frame(mock2)
colnames(mock2)[2] <- c("observed")

library(ggpubr)

l <- ggscatter(mock2, x = "expected", y = "observed", color = "Species", size = 5, alpha = 1,
               add = "reg.line",  add.params = list(color = "black", fill = "lightgray"),
               conf.int = TRUE) + stat_cor(aes(label = paste(..rr.label.., sep = "~`,`~")),
                                           label.x = 1) +
  geom_point(pch=21, size = 5)+
  theme(text = element_text(size = 12),
        legend.position = "right",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  scale_color_viridis(option = "rocket", discrete = T) +
  ylab("Observed0") + xlab("Esperado") +
  geom_hline(yintercept = 0.10, linetype="dotted") +
  ylim(0, 0.15) + labs(color = "Especie") +
  ggtitle("match")
l
############################################################################################################################################

list.files()

c1 <- as.data.frame(read.table("COUNTS_polished.tsv", sep = "\t", header = T))
t1 <- as.data.frame(read.table("TAXCLA_polished.tsv", sep = "\t", header = T))
name <- as.data.frame(colnames(c1))
df <- as.data.frame(str_split_fixed(name$`colnames(c1)`, "_", 3))
df1 <- as.data.frame(cbind(name, df$V2))
colnames(df1) <- c("file", "id")
head(df1)

meta1 <- df1
meta1$round <- c('polish')
row.names(meta1) <- meta1$file
head(meta1)
dim(meta1)


c1 <- as.matrix(read.table("COUNTS_polished.tsv", sep = "\t", header = T))
t1 <- as.matrix(read.table("TAXCLA_polished.tsv", sep = "\t", header = T))

OTU1 = otu_table(c1, taxa_are_rows = TRUE)
TAX1 = tax_table(t1)
samples = sample_data(meta1)
OTU1
TAX1

ps1 = phyloseq(OTU1, TAX1, samples)
ps1

ps <- tax_glom(ps1, taxrank = "Genus")
ps <- psmelt(ps1)

pp <- ps[which(ps$id == "DRR225046"),]

top20 <- names(sort(taxa_sums(ps1), decreasing=TRUE))[1:10000]
ps.top20 <- transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

graph1 <- psmelt(ps.top20)
head(graph1)
tab2 <- graph1
tab2$step <- c("polish")

tab2[grepl("Klebsiella", tab2$Species),]
tab2[grepl("Acinetobacter", tab2$Species),]
tab2[grepl("Enterococcus", tab2$Species),]
tab2[grepl("Enterobacter", tab2$Species),]
tab2[grepl("Pseudomonas", tab2$Species),]
tab2[grepl("Staphylococcus", tab2$Species),]


write.table(tab2, "polished_abundance_tab.tsv", sep = "\t", row.names = F, quote = F)

mock <- graph1[which(graph1$id == "DRR225046"),]

unique(mock$Species)

target <- c("Cereibacter sphaeroides", "Streptococcus mutans", "Enterococcus faecalis", "Bacillus cereus", "Clostridium beijerinckii",
            "Bifidobacterium adolescentis", "Escherichia coli", "Lactobacillus gasseri", "Staphylococcus epidermidis", "Streptococcus mutans", "Deinococcus radiodurans")

mock4 <- filter(mock, Species %in% target)
mock4 <- mock4 %>% group_by(Species) %>%
  summarise(sum = sum(Abundance))
max(mock4$sum)
min(mock4$sum)
unique(mock4$Species)

mock4$expected <- c("0.10")
mock4 <- as.data.frame(mock4)
colnames(mock4)[2] <- c("observed")

ll <- ggscatter(mock4, x = "expected", y = "observed", color = "Species", size = 5, alpha = 1,
                add = "reg.line",  add.params = list(color = "black", fill = "lightgray"),
                conf.int = TRUE) + stat_cor(aes(label = paste(..rr.label.., sep = "~`,`~")),
                                            label.x = 1) +
  geom_point(pch=21, size = 5)+
  theme(text = element_text(size = 12),
        legend.position = "right",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  scale_color_viridis(option = "rocket", discrete = T) +
  ylab("Observed0") + xlab("Esperado") + 
  geom_hline(yintercept = 0.10, linetype="dotted") + 
  ylim(0, 0.15) + labs(color = "Especie") +
  ggtitle("polish")

c <- ggarrange(l, ll, ncol = 2, common.legend = T, legend = "right")

png("/mnt/raid2tb/bubble/tesis/capitulo2/Figures/compare_mock_matsuo_porefile.png", res = 600, height = 20, width = 25, units = "cm")
c
dev.off()

