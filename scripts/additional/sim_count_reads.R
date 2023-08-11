library(phyloseq)
library(ggplot2)
library(stringr)
library(dplyr)
library(Metrics)
library(RColorBrewer)
library(viridis)

########################################
######Output pre-processing#############
########################################

setwd("/mnt/raid2tb/bubble/tesis/capitulo2/Fig3_S1")

setwd("./results_silva_zymo_noyacrd_simdefault/")
list.files()

############
##porefile##
############

#generate metadata from file names
c1 <- as.data.frame(read.table("COUNTS.tsv", sep = "\t", header = T))
t1 <- as.data.frame(read.table("TAXCLA.tsv", sep = "\t", header = T))
name <- as.data.frame(colnames(c1))
df <- as.data.frame(str_split_fixed(name$`colnames(c1)`, "_", 3))
df1 <- as.data.frame(cbind(name, df$V1))
colnames(df1) <- c("file", "id")
head(df1)

meta1 <- df1
meta1$round <- c('match')
row.names(meta1) <- meta1$file
head(meta1)
dim(meta1)

#load count table and taxonomy to a phyloseq object
c1 <- as.matrix(read.table("COUNTS.tsv", sep = "\t", header = T))
t1 <- as.matrix(read.table("TAXCLA.tsv", sep = "\t", header = T))

OTU1 = otu_table(c1, taxa_are_rows = TRUE)
TAX1 = tax_table(t1)
samples = sample_data(meta1)
OTU1
TAX1

ps1 = phyloseq(OTU1, TAX1, samples) #phyloseq object
ps1

rel <- psmelt(ps1) #get the read count per taxa

sp <- tax_glom(ps1, taxrank = "Species")
sp_tab <- psmelt(sp)
sp_tab$Rank <- c("Species")
mean(sp_tab$Abundance)

a <- sp_tab %>% group_by(Sample, Rank) %>%
  summarise(sum = sum(Abundance))

gn <- tax_glom(ps1, taxrank = "Genus")
gn_tab <- psmelt(gn)
gn_tab$Rank <- c("Genus")

mean(gn_tab$Abundance)

b <- gn_tab %>% group_by(Sample, Rank) %>%
  summarise(sum = sum(Abundance))


fm <- tax_glom(ps1, taxrank = "Family")
fm_tab <- psmelt(fm)
fm_tab$Rank <- c("Family")

mean(fm_tab$Abundance)

c <- fm_tab %>% group_by(Sample, Rank) %>%
  summarise(sum = sum(Abundance))


or <- tax_glom(ps1, taxrank = "Order")
or_tab <- psmelt(or)
or_tab$Rank <- c("Order")

mean(or_tab$Abundance)

d <- or_tab %>% group_by(Sample, Rank) %>%
  summarise(sum = sum(Abundance))


cl <- tax_glom(ps1, taxrank = "Class")
cl_tab <- psmelt(cl)
cl_tab$Rank <- c("Class")

mean(cl_tab$Abundance)


e <- cl_tab %>% group_by(Sample, Rank) %>%
  summarise(sum = sum(Abundance))


py <- tax_glom(ps1, taxrank = "Phylum")
py_tab <- psmelt(py)
py_tab$Rank <- c("Phylum")

mean(py_tab$Abundance)

f <- py_tab %>% group_by(Sample, Rank) %>%
  summarise(sum = sum(Abundance))


k <- tax_glom(ps1, taxrank = "Kingdom")
k_tab <- psmelt(k)
k_tab$Rank <- c("Kingdom")

mean(k_tab$Abundance)

g <- k_tab %>% group_by(Sample, Rank) %>%
  summarise(sum = sum(Abundance))


tab <- as.data.frame(rbind(a, b, c, d, e, f, g))
tab$Rank <- gsub("Kingdom", "Reino", tab$Rank)
tab$Rank <- gsub("Phylum", "Filo", tab$Rank)
tab$Rank <- gsub("Class", "Clase", tab$Rank)
tab$Rank <- gsub("Order", "Orden", tab$Rank)
tab$Rank <- gsub("Family", "Familia", tab$Rank)
tab$Rank <- gsub("Genus", "Género", tab$Rank)
tab$Rank <- gsub("Species", "Especie", tab$Rank)


level_order <- c("Reino", "Filo", "Clase", "Orden", "Familia", "Género", "Especie") 

p <- ggplot(tab, aes(x=Rank, y=sum)) + 
  geom_boxplot() +
  scale_x_discrete(limits = level_order) +
  scale_y_continuous(limits = c(25000, 40000)) +
  ylab("Número de lecturas clasificadas") +
  xlab("Grupo taxonómico") +
  theme_minimal()
p
