library(phyloseq)
library(ggplot2)
library(stringr)
library(dplyr)
library(Metrics)
library(viridis)
library(dplyr)
library(Metrics)
library(ggpubr)
library(data.table)

setwd("/mnt/raid2tb/bubble/tesis/capitulo2/Fig5_local/")


tab <- read.table("table_fig5.tsv", sep = "\t", header = T)


p <- ggplot(tab, aes(x=method_f, y=observed)) +
  geom_boxplot(alpha = 0) +
  theme(legend.position="bottom", 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        text = element_text(size = 40),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.text.x = element_text(size = 32),
        axis.text.y = element_text(size = 32),
        strip.text = element_text(size = 32),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  ylab("Abundance (%)") + 
  labs(fill="Source") +
  geom_signif(comparisons = list(c("porefile", "EMU")), map_signif_level = T, y_position = 25, size = 1, textsize = 15) +
  geom_signif(comparisons = list(c("porefile", "wf-metagenomics")), map_signif_level = T, y_position = 28, size = 1, textsize = 15) +
  geom_signif(comparisons = list(c("EMU", "wf-metagenomics")), map_signif_level = T, y_position = 27, size = 1, textsize = 15) +
  scale_fill_manual(values = c("#5c8001", "#f77f00")) + 
  scale_y_continuous(limits = c(0, 35)) +
  geom_hline(yintercept = 8.3, linetype = "dotted") +
  geom_hline(yintercept = 0.99, linetype = "dotted") +
  geom_jitter(size = 10, pch=21, aes(fill = host)) 

p

unique(tab$host)
tab$host <- gsub("human", "humano", tab$host)
tab$host <- gsub("environment", "ambiental", tab$host)

ls <- ggscatter(tab, x = "expected", y = "observed", color = "host", size = 10, alpha = 0.8,
                palette = c("#5c8001", "#f77f00"),
                add = "reg.line",  add.params = list(color = "black", fill = "lightgray"),
                conf.int = TRUE) + stat_cor(size =12, aes(label = paste(..rr.label.., sep = "~`,`~")),
                                            label.x = 1) +
  geom_point(pch=21, size = 10)+
  theme(text = element_text(size = 40),
        legend.position = "bottom",
        axis.text.x = element_blank(),
        axis.text.y = element_text(size = 24),
        strip.text = element_text(size = 32),
        strip.background = element_rect(color="black", fill="white", size=1, linetype="solid")) +
  labs(color="Orígen") +
  ylab("Observedo") + xlab("Esperado") + 
  scale_x_continuous(breaks=seq(0,30,by=10)) +
  geom_hline(yintercept = c(0.99, 8.3), linetype = "dotted") +
  facet_grid(~method_f, scales = "free_x") 
ls

###################################################################################################
###################################################################################################


head(tab)
tab$dis <- tab$expected - tab$observed
tab[is.na(tab)] <- 0

ptab <- select(tab, Species, expected, host)
ptab$method <- c("Expected")
colnames(ptab)[2] <- c("abundance")
head(ptab) 

ctab <- select(tab, Species, observed, method, host)
colnames(ctab)[2] <- c("abundance")
head(ctab)
c <- as.data.frame(rbind(ptab, ctab))
head(c)
c <- c %>% filter(abundance > 0)

c$method_f <- factor(c$method, levels=c("Expected", "porefile", "EMU", "wf-metagenomics"))
c$method <- NULL

wide <- reshape(c, idvar = c("Species", "host"), timevar = "method_f", direction = "wide")

s2 <- as.data.frame(wide)
colnames(s2) <- c("Especie", "Orígen", "Esperado", "porefile", "EMU", "wf-metagenomics")
head(s2)

write.table(s2, "tablaS2_c2.tsv", sep = "\t", row.names = F, quote = F)

#####################################################################################################
#####################################################################################################
l <- melt(setDT(wide), id.vars = c("Species", "host"), variable.name = "method")
ll <- l
l <- na.omit(l)

tab <-as.data.frame(table(l$method))
tab$pct <- tab$Freq/112*100
tab$type <- c("Total")

h <- l[which(l$host == "humano"),]
tabh <-as.data.frame(table(h$method))
tabh$pct <- tabh$Freq/100*100
tabh$type <- c("Humano")

e <- l[which(l$host == "ambiental"),]
tabe <-as.data.frame(table(e$method))
tabe$pct <- tabe$Freq/12*100
tabe$type <- c("Ambiental")


tabt <- as.data.frame(rbind(tab, tabh, tabe))
tabt <- tabt[-which(tabt$Var1 == "abundance.Expected"),]
tabt$Var1 <- gsub("abundance.", "", tabt$Var1)
tabt$Freq <- NULL
colnames(tabt) <- c("Metodo", "pct", "Tipo")
tabt$pct <- round(tabt$pct, digits = 0)

f <- ggplot(data=tabt, aes(x=Metodo, y=pct, fill = Metodo)) +
  geom_bar(stat="identity")+
  geom_text(aes(label=pct), size = 8, vjust=-0.3, size=3.5)+
  facet_grid(~Tipo, scales = "free_y") +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.y = element_text(size = 32),
        strip.text = element_text(size = 32),
        text = element_text(size = 40),
        strip.background = element_rect(color="black", fill="white", size=.1, linetype="solid")) +
  scale_fill_manual(values = c("porefile" = "#0d3b66", "EMU" = "#f95738", "wf-metagenomics" = "#fca311")) +
  ylab("Clasificación correcta (%)") +
  labs(fill = "Método")
f


ll <- ll %>% filter(if_any(everything(), is.na))
table(ll$method)

pp <- ll[which(ll$method == "abundance.porefile"),]
e <- ll[which(ll$method == "abundance.EMU"),]
w <- ll[which(ll$method == "abundance.wf-metagenomics"),]

intersect(intersect(pp$Species, e$Species), w$Species)


##############################################################################################
##############################################################################################

set.seed(227)
w <- wide[sample(nrow(wide), 20), ]

library(data.table)
long <- melt(setDT(w), id.vars = c("Species", "host"), variable.name = "method")
colnames(long) <- c("Species", "host", "method", "value")
long$method <- gsub("abundance.", "", long$method)
long$method <- gsub("Expected", "Esperado", long$method)

long$method_f <- factor(long$method, levels=c("Esperado", "porefile", "EMU", "wf-metagenomics"))
long <- long %>% filter(value > 0)
long <- as.data.frame(long)

library(RColorBrewer)
coul <- brewer.pal(8, "Spectral")
coul <- colorRampPalette(coul)(10)

head(long)
long$value <- as.numeric(long$value)


q <- ggplot(long, aes(x=method_f, y=Species, fill=value))+
  geom_tile(colour="white", size=0.5) +
  scale_fill_gradient2(low = "#f6f4d2", mid = "#f6bd60", high = "#e63946") +
  theme(legend.position="bottom", 
        text = element_text(size = 32),
        legend.text=element_text(size=12),
        legend.title=element_text(size=20),
        axis.line = element_line(colour = "black"), 
        axis.text.y = element_text(size = 32, face = "italic"),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_text(size = 32, face = "bold"),
        axis.text.x = element_text(size = 32, angle = 30, hjust = 1, vjust = 1),
        strip.background = element_rect(color="black", fill="white", size=.1, linetype="solid")) +
  labs(fill='Abundancia (%)') +
  ylab("Especie")
q

g <- ggarrange(ls, f, nrow = 2, legend = "bottom", labels = c("", "C"), 
               font.label = list(size = 28, color = "black"))
g

fig5 <- ggarrange(q, g, ncol = 2, heights = c(0.4, 0.6), widths = c(0.4, 0.6), labels = c("A", "B"), 
                  font.label = list(size = 28, color = "black"))

png("/mnt/raid2tb/bubble/tesis/capitulo2/Figuras/Figura5.png", res = 600, height = 50, width = 65, units = "cm")
fig5
dev.off()

