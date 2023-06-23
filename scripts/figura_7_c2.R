library(phyloseq)
library(ggplot2)
library(stringr)
library(dplyr)
library(Metrics)
library(viridis)

#esto esta aca /mnt/cive/csalazar/projects/3.porfile/real_data_public/matsuo/nanopore/filtered/samples


setwd("/mnt/raid2tb/bubble/tesis/capitulo2/Fig6/results/")
library(RColorBrewer)
coul <- brewer.pal(11, "Spectral")
coul <- colorRampPalette(coul)(50)
list.files()

c1 <- as.data.frame(read.table("COUNTS.tsv", sep = "\t", header = T))
t1 <- as.data.frame(read.table("TAXCLA.tsv", sep = "\t", header = T))
name <- as.data.frame(colnames(c1))
df <- as.data.frame(str_split_fixed(name$`colnames(c1)`, "_", 3))
df1 <- as.data.frame(cbind(name, df$V2))
colnames(df1) <- c("file", "id")
head(df1)

meta1 <- df1
meta1$round <- c('match')
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

top20 <- names(sort(taxa_sums(ps1), decreasing=TRUE))[1:50]
ps.top20 <- transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

graph1 <- psmelt(ps.top20)

#graph1 <- graph1[-which(graph1$id == "DRR225044"),]
graph1 <- graph1[-which(graph1$id == "DRR225046"),]

head(graph1)
tab1 <- graph1
tab1$step <- c("match")

graph1 <- graph1 %>% filter(Abundance > 0.0039)

p <- ggplot(graph1, aes(x = id, y = Abundance, fill = Species))
g1 <- p +  geom_bar(aes(), stat="identity", position="fill") + 
  scale_fill_manual(values = coul) +
  theme(legend.position="bottom", axis.line = element_line(colour = "black"), 
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
g1

#############################################################################################################################

setwd("/mnt/raid2tb/bubble/tesis/capitulo2/Fig6/results/")

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

top20 <- names(sort(taxa_sums(ps1), decreasing=TRUE))[1:10000]
ps.top20 <- transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

graph2 <- psmelt(ps.top20)
unique(graph2$id)

write.table(graph2, "polished_abundance_tab.tsv", sep = "\t", row.names = F, quote = F)


top20 <- names(sort(taxa_sums(ps1), decreasing=TRUE))[1:50]
ps.top20 <- transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

graph2 <- psmelt(ps.top20)
unique(graph2$id)
graph2 <- graph2[-which(graph2$id == "DRR225046"),]

p <- ggplot(graph2, aes(x = id, y = Abundance, fill = Species))
g2 <- p +  geom_bar(aes(), stat="identity", position="fill") + 
  scale_fill_manual(values = coul) +
  theme(legend.position="bottom", axis.line = element_line(colour = "black"), 
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
g2

tab2 <- graph2
tab2$step <- c("polish")


###############################################################################################################################################
top20 <- names(sort(taxa_sums(ps1), decreasing=TRUE))[1:10000]
ps.top20 <- transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

graph2 <- psmelt(ps.top20)
mock <- graph2[which(graph2$id == "DRR225046"),]
unique(mock$Species)

target <- c("Cereibacter sphaeroides", "Streptococcus mutans", "Enterococcus faecalis", "Bacillus cereus", "Clostridium beijerinckii",
            "Bifidobacterium adolescentis", "Escherichia coli", "Lactobacillus gasseri", "Staphylococcus epidermidis", "Streptococcus mutans", "Deinococcus radiodurans")

mock <- filter(graph2, Species %in% target)
mock <- mock %>% group_by(Species) %>%
  summarise(sum = sum(Abundance))
max(mock$sum)
min(mock$sum) # 0.03952569 es el límite de detección de la mock

##########################################################################################################################################################

tab <- as.data.frame(rbind(tab1, tab2))

tab <- tab %>%
  mutate(id2 = case_when(
    endsWith(id, "DRR225048") ~ "G1",
    endsWith(id, "DRR225051") ~ "G2",
    endsWith(id, "DRR225054") ~ "G3",
    endsWith(id, "DRR225057") ~ "G4",
    endsWith(id, "DRR225060") ~ "G5",
    endsWith(id, "DRR225063") ~ "G6"))

col <- c("#abaa30",
         "#df95ff",
         "#4bf548",
         "#e9a4ff",
         "#8de200",
         "#9a9ef1",
         "#b7e400",
         "#76a9ff",
         "#aee930",
         "#ff9bdd",
         "#62c200",
         "#e786bf",
         "#00c74d",
         "#ff82a1",
         "#5df18c",
         "#e8beff",
         "#a4c900",
         "#01b1fc",
         "#cbd100",
         "#a1baff",
         "#80b50b",
         "#81a6e7",
         "#ade856",
         "#ffa6cb",
         "#0ebe5a",
         "#f7847b",
         "#00e69c",
         "#ff9da1",
         "#01d07d",
         "#ffa087",
         "#02d3a1",
         "#ffb348",
         "#50d3ff",
         "#ffc03b",
         "#62e3ff",
         "#d59b16",
         "#01e6e7",
         "#ffb15f",
         "#01bccf",
         "#a3af00",
         "#56b0d5",
         "#94b10e",
         "#02c9b2",
         "#e29358",
         "#3bba83",
         "#f1d561",
         "#7ccfbb",
         "#d7de50",
         "#94e8b8",
         "#bfa345",
         "#99ea7e",
         "#f3d382",
         "#76b46f",
         "#dadb83",
         "#93b041",
         "#c7dea6",
         "#c3e27d",
         "#bebf86",
         "#bee292",
         "#a5aa65")

p <- ggplot(tab, aes(x = step, y = Abundance, fill = Species))
g3 <- p +  geom_bar(aes(), stat="identity", position="fill") + 
  scale_fill_manual(values = coul) +
  theme(legend.position="right", axis.line = element_line(colour = "black"), 
        text = element_text(size = 16), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 10, angle = 70, vjust = 0.8, hjust = 1),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  guides(fill=guide_legend(ncol = 2)) + ggtitle("") +
  facet_grid(~id2, scales = "free_x") +
  xlab("") + ylab("Abundancia relativa") + labs(fill = "Especie") 
g3



png("/mnt/raid2tb/bubble/tesis/capitulo2/Figuras/Figura7.png", res = 600, height = 20, width = 30, units = "cm")
g3
dev.off()

