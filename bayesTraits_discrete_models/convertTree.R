rm(list=ls())
setwd("D:/Ophiostoma/comparative_genomics/2023_01_sordariomycetes_revision/06_bayestraits_biggenomes/models_pathogenicity")
dir()

library(ape)
library(phytools)


### Original tree
tree <- read.nexus("filter_combine_1000_563.phy_ROOT.nexus")
tree
tree$tip.label



for (i in 0:9) {
  print(i)
  fname = paste0("../subsets/getSubset_",i,".txt")
  df <- read.table(fname, sep='\t', header=FALSE)
  species <- df$V1
  tree2 <- keep.tip(tree,tip=species)
  tree3 <- root(tree2, outgroup = c("ASARC","GLOZO","CLONG","HHEPA",
                                    "LHYAL","BCINE","SSCLE","PDEST","SKOCH",
                                    "TGUIA","XHEVE"),resolve.root = TRUE)
  write.nexus(tree3, file=paste0("./subset_",i,"/filter_combine_1000_SUBSET.phy_ROOT.nexus"),translate=TRUE)
}



### Writing tree
#print(paste0("filter_combine_1000_",length(good_samples),".phy_ROOT.nexus"))
#write.nexus(tree3, file=paste0("filter_combine_1000_",length(good_samples),".phy_ROOT.nexus"),translate=TRUE)
