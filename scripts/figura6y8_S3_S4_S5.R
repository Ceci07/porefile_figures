library(fs) # File System Operations Based on 'libuv'
library(magrittr)
library(vroom)
library(data.table) # Extension of `data.frame`
library(formattable) # Create 'Formattable' Data Structures
library(tidyverse)  # attaches purrr and readr
library(tidyselect)
library(writexl)
library(ggalluvial)
library(ggpubr)
library(grid)
library(gridExtra)
library(viridis)

setwd("/mnt/raid2tb/bubble/tesis/capitulo2/Fig6/emu_new/all/")
list.files()
# read file path
all_paths <-
  list.files(path = "./",
             pattern = "*.tsv",
             full.names = TRUE)
# read file content
all_content <-
  all_paths %>%
  lapply(read.table,
         header = TRUE,
         sep = "\t",
         encoding = "UTF-8")


# read file name
all_filenames <- all_paths %>%
  basename() %>%
  as.list()

# combine file content list and file name list
all_lists <- mapply(c, all_content, all_filenames, SIMPLIFY = FALSE)

# unlist all lists and change column name
all_result <- rbindlist(all_lists, fill = T)
# change column name
names(all_result)[14] <- "File.Path"
all_result$File.Path <- gsub("_rel-abundance.tsv", "", all_result$File.Path)
head(all_result)

emu <- select(all_result, File.Path, abundance, species, genus)
emu$abundance <- gsub("abundance", "", emu$abundance)
colnames(emu) <- c("id", "Abundance", "Species", "Genus")
emu$Abundance <- as.numeric(emu$Abundance)
emu$id <- gsub("filt_", "", emu$id)
emu <- as.data.frame(emu)
emu <- emu[complete.cases(emu),]
emu <- select(emu, id, Abundance, Species, Genus)
head(emu)

emuu <- emu
head(emuu)
emuu <- select(emuu, id, Abundance, Species, Genus)
#emuu <- emuu[which(emuu$id == "DRR225063"),]

emuu <- emuu %>% group_by(id, Species) %>%
  summarise(sum  = sum(Abundance))
emuu <- emuu %>% filter(sum > 0.01)
emuu$tool <- c("EMU")

library(dplyr)
tab1 <- read.table("/mnt/raid2tb/bubble/tesis/capitulo2/Fig6/results/polished_abundance_tab.tsv", 
                  sep = "\t", header = T)
head(tab1)
tabb <- tab1
tabb <- select(tabb, id, Abundance, Species, Genus)
tabb <- tabb %>% 
  group_by(id, Species) %>%
  summarise(sum  = sum(Abundance))
tabb$tool <- c("porefile")

tab <- as.data.frame(rbind(emuu, tabb))
tab <- tab %>% group_by(id, Species, tool) %>%
  summarise(sum  = sum(sum))

tab <- tab %>% filter(sum > 0.01)
# tab <- tab[-which(tab$id == "DRR225044"),]
tab <- tab[-which(tab$id == "DRR225046"),]
tab <- tab %>%
  mutate(id2 = case_when(
    endsWith(id, "DRR225048") ~ "G1",
    endsWith(id, "DRR225051") ~ "G2",
    endsWith(id, "DRR225054") ~ "G3",
    endsWith(id, "DRR225057") ~ "G4",
    endsWith(id, "DRR225060") ~ "G5",
    endsWith(id, "DRR225063") ~ "G6"))

coul <- brewer.pal(11, "Spectral")
coul <- colorRampPalette(coul)(70)

s <- ggplot(data = tab,
             aes(x = tool, y = sum, alluvium = Species, 
                 fill = Species, stratum = Species)) +
  geom_alluvium() +
  geom_stratum() +
  theme(panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        panel.background = element_blank(), 
        legend.position = "bottom",
        legend.text=element_text(size = 12),
        legend.title=element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 14),
        strip.background = element_rect(color="black", fill="white", size=0.5)) +
  scale_fill_manual(values = coul) +
  guides(fill=guide_legend(ncol=5, nrow = 15)) +
  facet_wrap(~id2, scales = "free_x", nrow = 2, ncol = 3)
s

ss <- annotate_figure(s, left = text_grob("Abundancia relativa", rot = 90, size =20))

png("/mnt/raid2tb/bubble/tesis/capitulo2/Figuras/FiguraS5.png", res = 600, height = 40, width = 45, units = 'cm')
ss
dev.off()

#################################################################################
emuu <- emu
head(emuu)
emuu <- select(emuu, id, Abundance, Species, Genus)
#emuu <- emuu[which(emuu$id == "DRR225063"),]

emuu <- emuu %>% group_by(id, Genus) %>%
  summarise(sum  = sum(Abundance))
emuu <- emuu %>% filter(sum > 0.01)
emuu$tool <- c("EMU")

head(tab1)
tabb <- tab1
tabb <- select(tabb, id, Abundance, Species, Genus)
tabb <- tabb %>% group_by(id, Genus) %>%
  summarise(sum  = sum(Abundance))
tabb$tool <- c("porefile")

tab <- as.data.frame(rbind(emuu, tabb))
tab <- tab %>% group_by(id, Genus, tool) %>%
  summarise(sum  = sum(sum))

tab <- tab %>% filter(sum > 0.01)
# tab <- tab[-which(tab$id == "DRR225044"),]
tab <- tab[-which(tab$id == "DRR225046"),]
tab <- tab %>%
  mutate(id2 = case_when(
    endsWith(id, "DRR225048") ~ "G1",
    endsWith(id, "DRR225051") ~ "G2",
    endsWith(id, "DRR225054") ~ "G3",
    endsWith(id, "DRR225057") ~ "G4",
    endsWith(id, "DRR225060") ~ "G5",
    endsWith(id, "DRR225063") ~ "G6"))

coul <- brewer.pal(11, "Spectral")
coul <- colorRampPalette(coul)(40)

g <- ggplot(data = tab,
            aes(x = tool, y = sum, alluvium = Genus, 
                fill = Genus, stratum = Genus)) +
  geom_alluvium() +
  geom_stratum() +
  theme(panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank(), 
        panel.background = element_blank(), 
        legend.position = "bottom",
        legend.text=element_text(size = 12),
        legend.title=element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 14),
        strip.background = element_rect(color="black", fill="white", size=0.2)) +
  scale_fill_manual(values = coul) +
  guides(fill=guide_legend(ncol=8, nrow = 5)) +
  facet_wrap(~id2, scales = "free_x", nrow = 2, ncol = 3)
g

gg <- annotate_figure(g, left = text_grob("Abundancia relativa", rot = 90, size =20))


png("/mnt/raid2tb/bubble/tesis/capitulo2/Figuras/FiguraS4.png", res = 600, height = 20, width = 35, units = 'cm')
gg
dev.off()

########################################################################################################
#######################################################################
#############################

emuu <- emu
head(emuu)
emuu <- select(emuu, id, Abundance, Species, Genus)

emuu <- emuu %>% group_by(id, Species) %>%
  summarise(sum  = sum(Abundance))
emuu <- emuu %>% filter(sum > 0.01)
emuu$tool <- c("EMU")

head(tab1)
tabb <- tab1
tabb <- select(tabb, id, Abundance, Species, Genus)
tabb <- tabb %>% group_by(id, Species) %>%
  summarise(sum  = sum(Abundance))
tabb$tool <- c("porefile")

tab <- as.data.frame(rbind(emuu, tabb))
tab <- tab %>% group_by(id, Species, tool) %>%
  summarise(sum  = sum(sum))
tab$Species <- gsub("Cereibacter sphaeroides", "Rhodobacter sphaeroides", tab$Species)

tab <- tab %>% filter(sum > 0.037)
tab <- tab[which(tab$id == "DRR225046"),]
class(tab)

library(RColorBrewer)
coul <- brewer.pal(11, "Spectral")
coul <- colorRampPalette(coul)(10)

target <- c("Rhodobacter sphaeroides", "Streptococcus mutans", "Enterococcus faecalis", "Bacillus cereus", "Clostridium beijerinckii",
            "Bifidobacterium adolescentis", "Escherichia coli", "Lactobacillus gasseri", "Staphylococcus epidermidis", "Streptococcus mutans",
            "Deinococcus radiodurans")

tab2 <- filter(tab, Species %in% target)
min(tab2$sum)

tab3 <- anti_join(tab, tab2, by = "Species")
tab3$Species <- NULL
tab3$Species <- c("Other")
tab3 <- select(tab3, id, Species, tool, sum)
tab <- as.data.frame(rbind(tab2, tab3))
head(tab)
tab$sum <- as.numeric(tab$sum)

s <- ggplot(data = tab,
            aes(x = tool, y = sum, alluvium = Species, 
                fill = Species, stratum = Species)) +
  geom_alluvium(alpha = 0.3) +
  geom_stratum() +
  theme(panel.grid.major = element_blank(), 
        axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 10), 
        axis.ticks.y = element_blank(), 
        panel.background = element_blank(), 
        legend.position = "bottom",
        legend.text=element_text(size = 12),
        legend.title=element_text(size=14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        strip.text = element_text(size = 14),
        strip.background = element_rect(color="black", fill="white", size=0.5)) +
  scale_fill_manual(values = coul) + ylab("Abundancia relativa") +
  guides(fill=guide_legend(ncol=3)) + labs(fill = "Especie") 
s

ss <- annotate_figure(s, left = text_grob("Abundancia relativa", rot = 90, size =20))

png("/mnt/raid2tb/bubble/tesis/capitulo2/Figuras/Figura6.png", res = 600, height = 20, width = 35, units = "cm")
ss
dev.off()



##################################################################################
####################################333
setwd("/home/csalazar/Documentos/tesis/capitulo2/parte1/matsuo/emu_new/all/")
list.files()
# read file path
all_paths <-
  list.files(path = "./",
             pattern = "*.tsv",
             full.names = TRUE)
# read file content
all_content <-
  all_paths %>%
  lapply(read.table,
         header = TRUE,
         sep = "\t",
         encoding = "UTF-8")


# read file name
all_filenames <- all_paths %>%
  basename() %>%
  as.list()

# combine file content list and file name list
all_lists <- mapply(c, all_content, all_filenames, SIMPLIFY = FALSE)

# unlist all lists and change column name
all_result <- rbindlist(all_lists, fill = T)
# change column name
names(all_result)[14] <- "File.Path"
all_result$File.Path <- gsub("_rel-abundance.tsv", "", all_result$File.Path)
head(all_result)

emu <- select(all_result, File.Path, abundance, family, genus)
emu$abundance <- gsub("abundance", "", emu$abundance)
colnames(emu) <- c("id", "Abundance", "Family", "Genus")
emu$Abundance <- as.numeric(emu$Abundance)
emu$id <- gsub("filt_", "", emu$id)
emu <- as.data.frame(emu)
emu <- emu[complete.cases(emu),]
head(emu)

emuu <- emu
head(emuu)
emuu <- select(emuu, id, Abundance, Genus)
emuu <- emuu %>% group_by(id, Genus) %>%
  summarise(sum  = sum(Abundance))
emuu$tool <- c("EMU")


tab1 <- read.table("/mnt/raid2tb/bubble/tesis/capitulo2/Fig6/results/polished_abundance_tab.tsv", 
                   sep = "\t", header = T)
head(tab1)
tabb <- tab1
tabb <- select(tabb, id, Abundance, Family, Genus)
tabb <- tabb %>% 
  group_by(id, Genus) %>%
  summarise(sum  = sum(Abundance))
tabb$tool <- c("porefile")


ill <- read.table("~/Documentos/porefile/real_dataset_matsuo/abundance_tables/illumina_abundance_matsuo.tsv", sep = "\t", header = T)
head(ill)
ill$id <- gsub(".fastq", "", ill$id)
ill <- select(ill, id, Family, Genus, Abundance)
ill <- ill %>% group_by(id, Genus) %>%
  summarise(sum  = sum(Abundance))
ill$tool <- c("DADA2")

tab <- as.data.frame(rbind(emuu, tabb, ill))
tab <- tab[-which(tab$id == "DRR225044"),]
tab <- tab[-which(tab$id == "DRR225046"),]

head(tab)
tab$Genus <- gsub("\\[|\\]", "", tab$Genus)
tab$Genus <- gsub("Ruminococcus torques group", "Ruminococcus", tab$Genus)
tab$Genus <- gsub("Ruminococcus gnavus group", "Ruminococcus", tab$Genus)
tab$Genus <- gsub("Clostridium innocuum group", "Clostridium", tab$Genus)
tab$Genus <- gsub("Clostridium sensu stricto 1", "Clostridium", tab$Genus)
tab$Genus <- gsub("Escherichia-Shigella", "Escherichia", tab$Genus)

tab <- tab %>% group_by(id, Genus, tool) %>%
  summarise(sum  = sum(sum))

tab <- tab %>% filter(sum > 0)

tab <- tab %>%
  mutate(id2 = case_when(
    endsWith(id, "DRR225048") ~ "G1",
    endsWith(id, "DRR225051") ~ "G2",
    endsWith(id, "DRR225054") ~ "G3",
    endsWith(id, "DRR225057") ~ "G4",
    endsWith(id, "DRR225060") ~ "G5",
    endsWith(id, "DRR225063") ~ "G6",
    endsWith(id, "DRR225050") ~ "G1",
    endsWith(id, "DRR225053") ~ "G2",
    endsWith(id, "DRR225056") ~ "G3",
    endsWith(id, "DRR225059") ~ "G4",
    endsWith(id, "DRR225062") ~ "G5",
    endsWith(id, "DRR225065") ~ "G6"))
unique(tab$Genus)
tab <- tab[complete.cases(tab), ]

q <- ggplot(tab, aes(x=tool, y=Genus, fill=sum))+
  geom_tile(colour="white", size=0.5) +
  scale_fill_viridis(option = "rocket", direction = -1) +
  labs(fill="AR") +
  theme(legend.position="right", axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(size = 5, face = "italic"),
        text = element_text(size = 16), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 12, angle = 60, hjust = 1, vjust = 1),
        strip.background = element_rect(color="black", fill="white", size=.1, linetype="solid")) +
  facet_grid(~id2, scales = "free_x") + ylab("Género")
q

png("/mnt/raid2tb/bubble/tesis/capitulo2/Figuras/FiguraS3.png", res = 600, height = 20, width = 20, units = "cm")
q
dev.off()

il <- tab[which(tab$tool == "DADA2"),]
il1 <- il[which(il$id2 == "G1"),]
il2 <- il[which(il$id2 == "G2"),]
il3 <- il[which(il$id2 == "G3"),]
il4 <- il[which(il$id2 == "G4"),]
il5 <- il[which(il$id2 == "G5"),]
il6 <- il[which(il$id2 == "G6"),]

emu <- tab[which(tab$tool == "EMU"),]
emu1 <- emu[which(emu$id2 == "G1"),]
emu2 <- emu[which(emu$id2 == "G2"),]
emu3 <- emu[which(emu$id2 == "G3"),]
emu4 <- emu[which(emu$id2 == "G4"),]
emu5 <- emu[which(emu$id2 == "G5"),]
emu6 <- emu[which(emu$id2 == "G6"),]

il_emu1 <- left_join(il1, emu1, by = "Genus")
il_emu1$status <- il_emu1$id2.x == il_emu1$id2.y

il_emu2 <- left_join(il2, emu2, by = "Genus")
il_emu2$status <- il_emu2$id2.x == il_emu2$id2.y

il_emu3 <- left_join(il3, emu3, by = "Genus")
il_emu3$status <- il_emu3$id2.x == il_emu3$id2.y

il_emu4 <- left_join(il4, emu4, by = "Genus")
il_emu4$status <- il_emu4$id2.x == il_emu4$id2.y

il_emu5 <- left_join(il5, emu5, by = "Genus")
il_emu5$status <- il_emu5$id2.x == il_emu5$id2.y

il_emu6 <- left_join(il6, emu6, by = "Genus")
il_emu6$status <- il_emu6$id2.x == il_emu6$id2.y

il_emu <- rbind(il_emu1, il_emu2, il_emu3, il_emu4, il_emu5, il_emu6)

il_emu$status <- as.character(il_emu$status)
il_emu$status <- il_emu$status %>% replace_na('not matched')
il_emu$status <- gsub("TRUE", "matched", il_emu$status)
table(il_emu$status)

porefile <- tab[which(tab$tool == "porefile"),]
porefile1 <- porefile[which(porefile$id2 == "G1"),]
porefile2 <- porefile[which(porefile$id2 == "G2"),]
porefile3 <- porefile[which(porefile$id2 == "G3"),]
porefile4 <- porefile[which(porefile$id2 == "G4"),]
porefile5 <- porefile[which(porefile$id2 == "G5"),]
porefile6 <- porefile[which(porefile$id2 == "G6"),]


il_porefile1 <- left_join(il1, porefile1, by = "Genus")
il_porefile1$status <- il_porefile1$id2.x == il_porefile1$id2.y

il_porefile2 <- left_join(il2, porefile2, by = "Genus")
il_porefile2$status <- il_porefile2$id2.x == il_porefile2$id2.y

il_porefile3 <- left_join(il3, porefile3, by = "Genus")
il_porefile3$status <- il_porefile3$id2.x == il_porefile3$id2.y

il_porefile4 <- left_join(il4, porefile4, by = "Genus")
il_porefile4$status <- il_porefile4$id2.x == il_porefile4$id2.y

il_porefile5 <- left_join(il5, porefile5, by = "Genus")
il_porefile5$status <- il_porefile5$id2.x == il_porefile5$id2.y

il_porefile6 <- left_join(il6, porefile6, by = "Genus")
il_porefile6$status <- il_porefile6$id2.x == il_porefile6$id2.y

il_porefile <- rbind(il_porefile1, il_porefile2, il_porefile3, il_porefile4, il_porefile5, il_porefile6)
il_porefile$status <- as.character(il_porefile$status)
il_porefile$status <- il_porefile$status %>% replace_na('not matched')
il_porefile$status <- gsub("TRUE", "matched", il_porefile$status)
table(il_porefile$status)


###############################################################################


# emu_porefile1 <- left_join(emu1, porefile1, by = "Genus")
# emu_porefile1$status <- emu_porefile1$id2.x == emu_porefile1$id2.y
# 
# emu_porefile2 <- left_join(emu2, porefile2, by = "Genus")
# emu_porefile2$status <- emu_porefile2$id2.x == emu_porefile2$id2.y
# 
# emu_porefile3 <- left_join(emu3, porefile3, by = "Genus")
# emu_porefile3$status <- emu_porefile3$id2.x == emu_porefile3$id2.y
# 
# emu_porefile4 <- left_join(emu4, porefile4, by = "Genus")
# emu_porefile4$status <- emu_porefile4$id2.x == emu_porefile4$id2.y
# 
# emu_porefile5 <- left_join(emu5, porefile5, by = "Genus")
# emu_porefile5$status <- emu_porefile5$id2.x == emu_porefile5$id2.y
# 
# emu_porefile6 <- left_join(emu6, porefile6, by = "Genus")
# emu_porefile6$status <- emu_porefile6$id2.x == emu_porefile6$id2.y
# 
# emu_porefile <- rbind(emu_porefile1, emu_porefile2, emu_porefile3, emu_porefile4, emu_porefile5, emu_porefile6)
# emu_porefile$status <- as.character(emu_porefile$status)
# emu_porefile$status <- emu_porefile$status %>% replace_na('not matched')
# emu_porefile$status <- gsub("TRUE", "matched", emu_porefile$status)
# table(emu_porefile$status)



#######################################################################################################################################

setwd("/home/csalazar/Documentos/tesis/capitulo2/parte1/matsuo/emu_new/all/")
list.files()
# read file path
all_paths <-
  list.files(path = "./",
             pattern = "*.tsv",
             full.names = TRUE)
# read file content
all_content <-
  all_paths %>%
  lapply(read.table,
         header = TRUE,
         sep = "\t",
         encoding = "UTF-8")


# read file name
all_filenames <- all_paths %>%
  basename() %>%
  as.list()

# combine file content list and file name list
all_lists <- mapply(c, all_content, all_filenames, SIMPLIFY = FALSE)

# unlist all lists and change column name
all_result <- rbindlist(all_lists, fill = T)
# change column name
names(all_result)[14] <- "File.Path"
all_result$File.Path <- gsub("_rel-abundance.tsv", "", all_result$File.Path)
head(all_result)

emu <- select(all_result, File.Path, abundance, family, genus)
emu$abundance <- gsub("abundance", "", emu$abundance)
colnames(emu) <- c("id", "Abundance", "Family", "Genus")
emu$Abundance <- as.numeric(emu$Abundance)
emu$id <- gsub("filt_", "", emu$id)
emu <- as.data.frame(emu)
emu <- emu[complete.cases(emu),]
head(emu)

emuu <- emu
head(emuu)
emuu <- select(emuu, id, Abundance, Genus)
emuu <- emuu %>% group_by(id, Genus) %>%
  summarise(sum  = sum(Abundance))
emuu$tool <- c("EMU")


tab1 <- read.table("/mnt/raid2tb/bubble/tesis/capitulo2/Fig6/results/polished_abundance_tab.tsv", 
                   sep = "\t", header = T)
head(tab1)
tabb <- tab1
tabb <- select(tabb, id, Abundance, Family, Genus)
tabb <- tabb %>% 
  group_by(id, Genus) %>%
  summarise(sum  = sum(Abundance))
tabb$tool <- c("porefile")


ill <- read.table("~/Documentos/porefile/real_dataset_matsuo/abundance_tables/illumina_abundance_matsuo.tsv", sep = "\t", header = T)
head(ill)
ill$id <- gsub(".fastq", "", ill$id)
ill <- select(ill, id, Family, Genus, Abundance)
ill <- ill %>% group_by(id, Genus) %>%
  summarise(sum  = sum(Abundance))
ill$tool <- c("DADA2")

tab <- as.data.frame(rbind(emuu, tabb, ill))
tab <- tab[-which(tab$id == "DRR225044"),]
tab <- tab[-which(tab$id == "DRR225046"),]

head(tab)
tab$Genus <- gsub("\\[|\\]", "", tab$Genus)
tab$Genus <- gsub("Ruminococcus torques group", "Ruminococcus", tab$Genus)
tab$Genus <- gsub("Ruminococcus gnavus group", "Ruminococcus", tab$Genus)
tab$Genus <- gsub("Clostridium innocuum group", "Clostridium", tab$Genus)
tab$Genus <- gsub("Clostridium sensu stricto 1", "Clostridium", tab$Genus)
tab$Genus <- gsub("Escherichia-Shigella", "Escherichia", tab$Genus)

tab <- tab %>% group_by(id, Genus, tool) %>%
  summarise(sum  = sum(sum))

tab2 <- tab %>% filter(sum > 0.005)

tab2 <- tab2 %>%
  mutate(id2 = case_when(
    endsWith(id, "DRR225048") ~ "G1",
    endsWith(id, "DRR225051") ~ "G2",
    endsWith(id, "DRR225054") ~ "G3",
    endsWith(id, "DRR225057") ~ "G4",
    endsWith(id, "DRR225060") ~ "G5",
    endsWith(id, "DRR225063") ~ "G6",
    endsWith(id, "DRR225050") ~ "G1",
    endsWith(id, "DRR225053") ~ "G2",
    endsWith(id, "DRR225056") ~ "G3",
    endsWith(id, "DRR225059") ~ "G4",
    endsWith(id, "DRR225062") ~ "G5",
    endsWith(id, "DRR225065") ~ "G6"))
unique(tab2$Genus)
tab2 <- tab2[complete.cases(tab2), ]
dim(as.data.frame(unique(tab2$Genus)))

q <- ggplot(tab2, aes(x=tool, y=Genus, fill=sum))+
  geom_tile(colour="white", size=0.5) +
  scale_fill_viridis(option = "rocket", direction = -1) +
  labs(fill="AR") +
  theme(legend.position="right", axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(size = 12, face = "italic"),
        text = element_text(size = 16), 
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 12, angle = 60, hjust = 1, vjust = 1),
        strip.background = element_rect(color="black", fill="white", size=.1, linetype="solid")) +
  facet_grid(~id2, scales = "free_x") + ylab("Género")
q

png("/mnt/raid2tb/bubble/tesis/capitulo2/Figuras/Figura8.png", res = 600, height = 20, width = 20, units = "cm")
q
dev.off()

