library(dplyr)
library(stringr)
library(ggplot2)
library(ggsignif)

setwd("./")

match <- read.table("matched_abundance_tab.tsv", sep = "\t", header = T)
polished <- read.table("polished_abundance_tab.tsv", sep = "\t", header = T)

head(match)
m <- match %>% group_by(id, Genus) %>%
  summarise(sum = sum(Abundance))
m$type <- c("match")

p <- polished %>% group_by(id, Genus) %>%
  summarise(sum = sum(Abundance))
p$type <- c("polish")

tab <- as.data.frame(rbind(m, p))
tab <- tab[order(tab$sum, decreasing = T),]

m2 <- m %>% filter(sum > 0.001)
p2 <- p %>% filter(sum > 0.001)

tab <- as.data.frame(rbind(m, p))

p <- ggplot(tab, aes(x=type, y=sum)) + 
  geom_boxplot() +
  facet_wrap(~Genus, nrow = 10, ncol = 9)+
  geom_signif(comparisons = list(c("match", "polish")), test= "wilcox.test", map_signif_level = T, y_position = c(0.7)) +
  theme_minimal()+
  ylab("Abundancia relativa") +
  geom_jitter() +
  xlab("")
p

png("genus_assess.png", res = 600, height = 40, width = 45, units = "cm")
p
dev.off()


#################################################

m <- match %>% group_by(id, Species) %>%
  summarise(sum = sum(Abundance))
m$type <- c("match")

p <- polished %>% group_by(id, Species) %>%
  summarise(sum = sum(Abundance))
p$type <- c("polish")

m2 <- m %>% filter(sum > 0.005)
p2 <- p %>% filter(sum > 0.005)

tab <- as.data.frame(rbind(m2, p2))

shapiro.test(tab$sum)

p <- ggplot(tab, aes(x=type, y=sum)) + 
  geom_boxplot() +
  facet_wrap(~Species, nrow = 10, ncol = 10)+
  geom_signif(comparisons = list(c("match", "polish")), test= "wilcox.test", map_signif_level = T, y_position = c(0.7)) +
  theme_minimal()+
  theme(text = element_text(size = 20)) +
  ylab("Abundancia relativa") +
  geom_jitter() +
  xlab("")
p

#################################################

m <- match %>% group_by(id, Family) %>%
  summarise(sum = sum(Abundance))
m$type <- c("match")

p <- polished %>% group_by(id, Family) %>%
  summarise(sum = sum(Abundance))
p$type <- c("polish")

m2 <- m %>% filter(sum > 0.001)
p2 <- p %>% filter(sum > 0.001)

tab <- as.data.frame(rbind(m2, p2))


shapiro.test(tab$sum)

p <- ggplot(tab, aes(x=type, y=sum)) + 
  geom_boxplot() +
  facet_wrap(~Family, nrow = 10, ncol = 9)+
  geom_signif(comparisons = list(c("match", "polish")), test= "wilcox.test", map_signif_level = T, y_position = c(0.7)) +
  theme_minimal()+
  ylab("Abundancia relativa") +
  geom_jitter() +
  xlab("")
p

#####################################################

m <- match %>% group_by(id,Phylum) %>%
  summarise(sum = sum(Abundance))
m$type <- c("match")

p <- polished %>% group_by(id, Phylum) %>%
  summarise(sum = sum(Abundance))
p$type <- c("polish")

tab <- as.data.frame(rbind(m, p))
tab <- tab[order(tab$sum, decreasing = T),]

m2 <- m %>% filter(sum > 0.001)
p2 <- p %>% filter(sum > 0.001)

tab <- as.data.frame(rbind(m2, p2))
# tab  <- tab  %>% drop_na()

p <- ggplot(tab, aes(x=type, y=sum)) + 
  geom_boxplot() +
  facet_wrap(~Phylum, nrow = 10, ncol = 9)+
  geom_signif(comparisons = list(c("match", "polish")), test= "wilcox.test", map_signif_level = T, y_position = c(0.7)) +
  theme_minimal()+
  ylab("Abundancia relativa") +
  geom_jitter() +
  xlab("")
p
