library(phyloseq)
library(ggplot2)
library(stringr)
library(dplyr)
library(RColorBrewer)
library(Metrics)
library(viridis)

setwd("./")

list.files()

c1 <- as.data.frame(read.table("COUNTS.tsv", sep = "\t", header = T))
t1 <- as.data.frame(read.table("TAXCLA.tsv", sep = "\t", header = T))
name <- as.data.frame(colnames(c1))
df <- as.data.frame(str_split_fixed(name$`colnames(c1)`, "_", 2))
dff <- as.data.frame(str_split_fixed(df$V2, "_", 2))
df1 <- as.data.frame(cbind(name, df$V2, df$V1))
colnames(df1) <- c("file", "id", "identity")
head(df1)
df1$identity <- gsub("X", "", df1$identity)

meta1 <- df1
meta1$round <- c('match')
row.names(meta1) <- meta1$file
head(meta1)
dim(meta1)
meta1 <- meta1[-which(meta1$identity == "80"),]

c1 <- as.matrix(read.table("COUNTS.tsv", sep = "\t", header = T))
t1 <- as.matrix(read.table("TAXCLA.tsv", sep = "\t", header = T))

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
tt <- psmelt(ps.top20)

target <- c("Listeria monocytogenes", "Pseudomonas aeruginosa", "Bacillus subtilis", "Escherichia coli",
            "Salmonella enterica","Limosilactobacillus fermentum", "Enterococcus faecalis", "Staphylococcus aureus")

int <- filter(tt, Species %in% target)
neg <- tt[!tt$Species %in% int$Species, ]
mean(neg$Abundance)
bs <- neg[grepl("Bacillus", neg$Species),]
unique(bs$Species)
min(int$Abundance)

top20 <- names(sort(taxa_sums(ps1), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

graph1 <- psmelt(ps.top20)
head(graph1)
tab1 <- graph1
tab1$step <- c("match")
unique(graph1$Species)


library(RColorBrewer)
coul <- brewer.pal(8, "Spectral")
coul <- colorRampPalette(coul)(20)

head(graph1)
p <- ggplot(graph1, aes(x = Sample, y = Abundance, fill = Species))
g <- p +  geom_bar(aes(), stat="identity", position="fill") + 
  scale_fill_manual(values = coul, name = "Especies") +
  facet_grid(~identity, scales = "free_x") +
  theme(legend.position="bottom", axis.line = element_line(colour = "black"), 
        text = element_text(size = 10), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  guides(fill=guide_legend(nrow=8, ncol = 8)) + ggtitle("") +
  xlab("") + ylab("Abundancia relativa") 
g
match <- g

#porefile read count
p_rc <- read.table("NanoPlots/summary.tsv", sep = "\t", header = T)
p <- psmelt(ps1)
t <- tax_glom(ps1, taxrank = "Species")
porefile <- psmelt(t)
porefile$method <- c("porefile")
head(porefile)

p <- porefile %>% group_by(id, identity) %>%
  summarise(sum = sum(Abundance))

p$FN <- 40000-p$sum
name <- as.data.frame(str_c(p$id, "-", p$identity))
colnames(name) <- c("name")
p <- as.data.frame(cbind(p, name))

target <- c("Listeria monocytogenes", "Pseudomonas aeruginosa", "Bacillus subtilis", "Escherichia coli",
            "Salmonella enterica","Limosilactobacillus fermentum", "Enterococcus faecalis", "Staphylococcus aureus")

int <- filter(porefile, Species %in% target)
TP <- int %>% group_by(id, identity) %>%
  summarise(sum = sum(Abundance))
colnames(TP) <- c("id", "ident", "TP")
name <- as.data.frame(str_c(TP$id, "-", TP$ident))
colnames(name) <- c("name")
TP <- as.data.frame(cbind(TP, name))
TP <- select(TP, TP, name)
p <- merge(p, TP, by = "name")

neg <- porefile[!porefile$Species %in% int$Species,]
FP <- neg %>% group_by(id, identity) %>%
  summarise(FP = sum(Abundance))
name <- as.data.frame(str_c(FP$id, "-", FP$identity))
colnames(name) <- c("name")
FP <- as.data.frame(cbind(FP, name))
FP <- select(FP, FP, name)
p <- merge(p, FP, by = "name")

p$precision <- p$TP/(p$TP + p$FP)*100
p$recall <- p$TP/(p$TP + p$FN)*100

mean(p$precision)
sd(p$precision)

mean(p$recall)
sd(p$recall)

p$method <- c("match")

percision <- p %>% group_by(identity) %>%
  summarise(mean = mean(precision), sd = sd(precision))
recall <- p %>% group_by(identity) %>%
  summarise(mean = mean(recall), sd = sd(recall))

pp <- select(p, name, FN, TP, FP)
ps <- pp[,-1]/40000*100
id <- as.data.frame(str_split_fixed(pp$name, "-", 2))
ps$name <- id$V2
pq1 <- ps %>% group_by(name) %>%
  summarise(meanFN = mean(FN), meanTP = mean(TP), meanFP = mean(FP))

pq2 <- ps %>% group_by(name) %>%
  summarise(sdFN = sd(FN), sdTP = sd(TP), sdFP = sd(FP))



##################################################################################
########################################################
#####################################

list.files()

c1 <- as.data.frame(read.table("COUNTS_polished.tsv", sep = "\t", header = T))
t1 <- as.data.frame(read.table("TAXCLA_polished.tsv", sep = "\t", header = T))
name <- as.data.frame(colnames(c1))
df <- as.data.frame(str_split_fixed(name$`colnames(c1)`, "_", 2))
dff <- as.data.frame(str_split_fixed(df$V2, "_", 2))
df1 <- as.data.frame(cbind(name, df$V2, df$V1))
colnames(df1) <- c("file", "id", "identity")
head(df1)
df1$identity <- gsub("X", "", df1$identity)

meta1 <- df1
meta1$round <- c('match')
row.names(meta1) <- meta1$file
head(meta1)
dim(meta1)
meta1 <- meta1[-which(meta1$identity == "80"),]

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
tt <- psmelt(ps.top20)

target <- c("Listeria monocytogenes", "Pseudomonas aeruginosa", "Bacillus subtilis", "Escherichia coli",
            "Salmonella enterica","Limosilactobacillus fermentum", "Enterococcus faecalis", "Staphylococcus aureus")

int <- filter(tt, Species %in% target)
neg <- tt[!tt$Species %in% int$Species, ]
mean(neg$Abundance)
bs <- neg[grepl("Bacillus", neg$Species),]
unique(bs$Species)
min(int$Abundance)

top20 <- names(sort(taxa_sums(ps1), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps1, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)

graph1 <- psmelt(ps.top20)
head(graph1)
tab1 <- graph1
tab1$step <- c("match")
unique(graph1$Species)


library(RColorBrewer)
coul <- brewer.pal(8, "Spectral")
coul <- colorRampPalette(coul)(20)

head(graph1)
p <- ggplot(graph1, aes(x = Sample, y = Abundance, fill = Species))
g <- p +  geom_bar(aes(), stat="identity", position="fill") + 
  scale_fill_manual(values = coul, name = "Especies") +
  facet_grid(~identity, scales = "free_x") +
  theme(legend.position="bottom", axis.line = element_line(colour = "black"), 
        text = element_text(size = 10), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_blank(),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  guides(fill=guide_legend(nrow=8, ncol = 8)) + ggtitle("") +
  xlab("") + ylab("Abundancia relativa") 
g
match <- g

#porefile read count
p_rc <- read.table("NanoPlots/summary.tsv", sep = "\t", header = T)
p <- psmelt(ps1)
t <- tax_glom(ps1, taxrank = "Species")
porefile <- psmelt(t)
porefile$method <- c("porefile")
head(porefile)

p <- porefile %>% group_by(id, identity) %>%
  summarise(sum = sum(Abundance))

p$FN <- 40000-p$sum
name <- as.data.frame(str_c(p$id, "-", p$identity))
colnames(name) <- c("name")
p <- as.data.frame(cbind(p, name))

target <- c("Listeria monocytogenes", "Pseudomonas aeruginosa", "Bacillus subtilis", "Escherichia coli",
            "Salmonella enterica","Limosilactobacillus fermentum", "Enterococcus faecalis", "Staphylococcus aureus")

int <- filter(porefile, Species %in% target)
TP <- int %>% group_by(id, identity) %>%
  summarise(sum = sum(Abundance))
colnames(TP) <- c("id", "ident", "TP")
name <- as.data.frame(str_c(TP$id, "-", TP$ident))
colnames(name) <- c("name")
TP <- as.data.frame(cbind(TP, name))
TP <- select(TP, TP, name)
p <- merge(p, TP, by = "name")

neg <- porefile[!porefile$Species %in% int$Species,]
FP <- neg %>% group_by(id, identity) %>%
  summarise(FP = sum(Abundance))
name <- as.data.frame(str_c(FP$id, "-", FP$identity))
colnames(name) <- c("name")
FP <- as.data.frame(cbind(FP, name))
FP <- select(FP, FP, name)
p <- merge(p, FP, by = "name")

p$precision <- p$TP/(p$TP + p$FP)*100
p$recall <- p$TP/(p$TP + p$FN)*100

mean(p$precision)
sd(p$precision)

mean(p$recall)
sd(p$recall)

p$method <- c("match")

percision <- p %>% group_by(identity) %>%
  summarise(mean = mean(precision), sd = sd(precision))
recall <- p %>% group_by(identity) %>%
  summarise(mean = mean(recall), sd = sd(recall))

pp <- select(p, name, FN, TP, FP)
ps <- pp[,-1]/40000*100
id <- as.data.frame(str_split_fixed(pp$name, "-", 2))
ps$name <- id$V2
pq1 <- ps %>% group_by(name) %>%
  summarise(meanFN = mean(FN), meanTP = mean(TP), meanFP = mean(FP))

pq2 <- ps %>% group_by(name) %>%
  summarise(sdFN = sd(FN), sdTP = sd(TP), sdFP = sd(FP))

