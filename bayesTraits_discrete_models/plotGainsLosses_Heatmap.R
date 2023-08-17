rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(RColorBrewer)

# Gains and losses
df <- read.csv('getTransitions.tsv', sep='\t', header=TRUE)
head(df)
unique(df$trait)

# Filtering out effectors
df <- df %>% filter(trait != 'effectors')

dg <- df %>% gather(key = "type", value = "fraction", c("gains", "no_changes","losses"))
head(dg)

dg %>% filter(trait == 'assembly_length')

# Model bayes factors
bf <- read.csv('getBayesFactor.out', sep='\t', header=TRUE)
head(bf)
unique(bf$trait)
dconvert <- data.frame("index" = c('00_01', '00_10', '01_11', '10_11'),
                       "rate" = c('mean_logBF_q12.q21', 'mean_logBF_q13.q31', 
                                         'mean_logBF_q24.q42', 'mean_logBF_q34.q43'))
dconvert

bf_long <- bf %>% dplyr::select(trait, mean_logBF_q12.q21, mean_logBF_q13.q31,
                     mean_logBF_q24.q42, mean_logBF_q34.q43, feature) %>% distinct() %>%
  gather(key = "rate", value = "BF", -trait, -feature)
bf_long <- merge(bf_long, dconvert, by = "rate", sort=FALSE)
bf_long$significance <- ifelse(bf_long$BF > 4.0, "1", "0")
head(bf_long)

# Merging

dg_merged <- merge(dg, bf_long, by = c("index", "trait"), sort=FALSE)
head(dg_merged)
unique(dg_merged$feature)

dg_merged %>% filter(trait == 'assembly_length')

#dg$significance <- ifelse((dg$index == "00_01" & dg$trait == "assembly_length"), "1", "0")
#head(dg)
dg_merged$type <- factor(dg_merged$type, levels = c("gains", "no_changes", "losses"), labels = c("Gain", "No change", "Loss"))
dg_merged$feature <- factor(dg_merged$feature, levels = c("Genome", "Genome w/o repeats", "Genes", "Repeats",
                                                          "Genes with introns", "Exon length", 
                                                          "Intergenic length", "Intron length", "Introns",
                                                          "GC", "Pseudo tRNA", " tRNA"))


p1 <- ggplot(dg_merged, aes(x = type, y = index, fill = fraction, alpha=significance)) +
  geom_tile(linewidth = 1.5,
            linetype = 1) +
  scale_alpha_manual(values = c(0.25, 1), guide=FALSE) +
  #scale_fill_distiller(direction=1) +
  scale_fill_gradient(low = "#A7C7E7", high = "#191970") + 
  geom_text(aes(label = round(fraction,2)), color = ifelse(dg_merged$fraction > 0.5, "white", "black"), size = 4) +
  coord_fixed() +
  labs(x = "Transition fraction", y = "Transition type") +
  facet_wrap(~feature, ncol=4) + 
  theme(panel.background = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size=14),
        legend.position = "none",
        axis.text = element_text(size=10),
        axis.title = element_text(size=14))



png("plotGainsLosses_Heatmap.png", w=1400, h = 1400, res=150)
p1
dev.off()
