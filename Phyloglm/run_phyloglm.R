rm(list=ls())

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
library(phylolm)
library(rr2)
library(phytools)
library(geiger)
library(adephylo)




### FUNCTIONS

# Rescaling tree
rescale_Tree <- function(tree,scale){
  tree$edge.length<-
    tree$edge.length/max(nodeHeights(tree)[,2])*scale
  return(tree)
}



### USING PHYLOGLM MODEL TO FIT GENOMIC TRAITS TO BINARY TRAIT: PATHOGENICITY (1/0)

for (sub in 0:9) {
  print(sub)
  
  # Getting dataset
  fname = paste0("../../06_bayestraits_biggenomes/subsets/getSubset_",sub,".csv")
  df <- read.table(fname, sep=',',header=TRUE)
  df$name <- df$code
  df$assembly_noreps <- df$assembly_length - df$repeats_length
  dff <- df
  traits <- c("assembly_length","genes","frac_repeats_dec","assembly_noreps","propGC",
              "N_introns","F_intron_genes","LEN_introns","LEN_exons","LEN_intergenic",
              "tRNA","pseudo_tRNA", "effectors")
  rownames(dff) <- dff$name
  dff$assembly_length_10M <- dff$assembly_length/10000000
  dff$genes_K <- dff$genes/1000
  dff$assembly_noreps_10M <- dff$assembly_noreps/10000000
  dff$LEN_introns_K <- dff$LEN_introns/1000
  dff$LEN_exons_100B <- dff$LEN_exons/100
  dff$LEN_intergenic_K <- dff$LEN_intergenic/1000
  dff$tRNA_100B <- dff$tRNA/100
  dff$pseudo_tRNA_100B <- dff$pseudo_tRNA/100
  dff$effectors_100B <- dff$effectors/100
  
  # Removing species with missing pathogenicity values
  unknown <- dff %>% filter(AnyPathogen == 'NAN') %>% pull(code) %>% as.vector()
  df2 <- dff %>% filter(!code %in% unknown)
  dim(df2)
  pathogens_binary <- as.numeric(as.character(df2$AnyPathogen))
  
  # Getting phylogeny & and selecting sampled species
  tree <- read.tree("tree_r8s_ultrametric_noAnchor.txt")
  tree2 <- keep.tip(tree,tip=rownames(df2))
  tree3 <- root(tree2, outgroup = c("ASARC","GLOZO","CLONG","HHEPA",
                                    "LHYAL","BCINE","SSCLE","PDEST","SKOCH",
                                    "TGUIA","XHEVE"),resolve.root = TRUE)
  treeS3 <- rescale_Tree(tree3, 1.0)
  
  
  # Running phyloglm for prepared traits
  traits <- c("assembly_length_10M","genes_K","frac_repeats_dec","assembly_noreps_10M","propGC",
              "N_introns","F_intron_genes","LEN_introns_K","LEN_exons_100B","LEN_intergenic_K",
              "tRNA_100B","pseudo_tRNA_100B","effectors_100B")
  method2 <- "logistic_MPLE"
  
  comboR <- c()

  for (i in 1:length(traits)) {
    trait <- traits[i]
    print(trait)
    f <- as.formula(paste("pathogens_binary",trait, sep=" ~ "))
    fit.phyloglm <- phylolm::phyloglm(f, data=df2, phy = treeS3, btol=30, boot=1000, method = method2)
    dt <- as.data.frame(summary(fit.phyloglm)$coefficients)
    rownames(dt) <- c()
    dt['trait'] <- trait
    dt['method'] <- method2
    dt['factor'] <- c("Intercept","trait")
    dt['logLik'] <- logLik.phyloglm(fit.phyloglm)$logLik
    
    if (exists('comboR')) {
      combo = bind_rows(comboR,dt)
      comboR = combo
    }
    else {
      comboR = dt
    }
  }
  
  head(comboR)
  comboR %>% filter(factor == 'trait')
  
  # Adjusting p-values
  comboM  <- comboR %>% filter(factor == 'trait') %>% arrange(trait)
  comboM$adjusted.p <- p.adjust(comboM$p.value)
  
  # Writing table
  write.table(comboM, file = paste0("run_phyloglm_results_",sub,".tab"), sep = "\t", col.names = TRUE, row.names = FALSE, quote=FALSE)
  
}












head(comboR)
comboR %>% filter(factor == 'trait')

### Calculating LRT
fit.1 <- phylolm::phyloglm(AnyPathogen ~ 1, data=df2, phy = tree2, btol=30, boot=1000, method = method1)
fit.2 <- phylolm::phyloglm(AnyPathogen ~ 1, data=df2, phy = tree2, btol=30, boot=1000, method = method2)
logs <- data.frame("method" = c(method1, method2),
                   "logLik0" = c(logLik.phyloglm(fit.1)$logLik, logLik.phyloglm(fit.2)$logLik))
comboM <- merge(comboR, logs, by = "method", sort=FALSE)
comboM['LR'] <- 2*(comboM$logLik - comboM$logLik0)
comboM['LRT_pval'] <- pchisq(comboM$LR, df = 1, lower.tail = FALSE)
head(comboM)


comboM  <- comboR %>% filter(factor == 'trait') %>% arrange(trait)

comboM$adjusted.p <- p.adjust(comboM$p.value)
# Writing table
write.table(comboM, file = "run_phyloglm_results.tab", sep = "\t", col.names = TRUE, row.names = FALSE, quote=FALSE)



