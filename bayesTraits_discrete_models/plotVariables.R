rm(list=ls())
setwd("D:/Ophiostoma/comparative_genomics/2023_01_sordariomycetes_revision/06_bayestraits_biggenomes/models_pathogenicity")
dir()

library(ggplot2)
library(dplyr)
library(tidyr)
library(cowplot)
#library(ggpubr)
library(tidyr)


for (i in 0:9) {
  print(i)
  fname = paste0("../subsets/getSubset_",i,".csv")
  df <- read.table(fname, sep=',', header=TRUE)
  ds1 <- df %>% select(filename,code,order,buscoCS,assembly_length,contigs,repeats_length,frac_repeats_dec,
                      genes,propGC,effectors,plant_pathogen,insect,outgroup,Plant.Saprobe,AnyPathogen,
                      N_introns,F_intron_genes,LEN_introns,N_exons,LEN_exons,LEN_intergenic,tRNA,pseudo_tRNA)
  ds1$assembly_noreps <- ds1$assembly_length - ds1$repeats_length
  
  ### Medians
  med_assembly_length <- median(ds1$assembly_length)
  med_genes <- median(ds1$genes)
  med_frac_repeats <- median(ds1$frac_repeats_dec)
  med_propGC <- median(ds1$propGC)
  med_assembly_noreps <- median(ds1$assembly_noreps)
  med_N_introns <- median(ds1$N_introns)
  med_F_intron_genes <- median(ds1$F_intron_genes)
  med_LEN_introns <- median(ds1$LEN_introns)
  med_LEN_exons <- median(ds1$LEN_exons)
  med_LEN_intergenic <- median(ds1$LEN_intergenic)
  med_tRNA <- median(ds1$tRNA)
  med_pseudo_tRNA <- median(ds1$pseudo_tRNA)
  med_effectors <- median(ds1$effectors)
  
  ds1$assembly_length_bin <- ifelse(ds1$assembly_length >= med_assembly_length,"1","0")
  ds1$genes_bin <- ifelse(ds1$genes >= med_genes,"1","0")
  ds1$frac_repeats_bin <- ifelse(ds1$frac_repeats_dec >= med_frac_repeats,"1","0")
  ds1$propGC_bin <- ifelse(ds1$propGC >= med_propGC,"1","0")
  ds1$assembly_noreps_bin <- ifelse(ds1$assembly_noreps >= med_assembly_noreps,"1","0")
  ds1$N_introns_bin <- ifelse(ds1$N_introns >= med_N_introns,"1","0")
  ds1$F_intron_genes_bin <- ifelse(ds1$F_intron_genes >= med_F_intron_genes,"1","0")
  ds1$LEN_introns_bin <- ifelse(ds1$LEN_introns >= med_LEN_introns,"1","0")
  ds1$LEN_exons_bin <- ifelse(ds1$LEN_exons >= med_LEN_exons,"1","0")
  ds1$LEN_intergenic_bin <- ifelse(ds1$LEN_intergenic >= med_LEN_intergenic,"1","0")
  ds1$tRNA_bin <- ifelse(ds1$tRNA >= med_tRNA,"1","0")
  ds1$pseudo_tRNA_bin <- ifelse(ds1$pseudo_tRNA >= med_pseudo_tRNA,"1","0")
  ds1$effectors_bin <- ifelse(ds1$effectors >= med_effectors,"1","0")
  
  ds1$patho0 <- ifelse(ds1$AnyPathogen == 1,"1","0")
  ds1$patho1 <- ifelse(ds1$AnyPathogen == "NAN","-",ds1$patho0)
  
  write.table(ds1[,c("code","genes_bin","patho1")],file=paste0("./subset_",i,"/bayestraits_genes.txt"),append=FALSE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(ds1[,c("code","assembly_length_bin","patho1")],file=paste0("./subset_",i,"/bayestraits_assembly_length.txt"),append=FALSE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(ds1[,c("code","frac_repeats_bin","patho1")],file=paste0("./subset_",i,"/bayestraits_frac_repeats.txt"),append=FALSE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(ds1[,c("code","propGC_bin","patho1")],file=paste0("./subset_",i,"/bayestraits_propGC.txt"),append=FALSE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(ds1[,c("code","assembly_noreps_bin","patho1")],file=paste0("./subset_",i,"/bayestraits_assembly_noreps.txt"),append=FALSE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(ds1[,c("code","N_introns_bin","patho1")],file=paste0("./subset_",i,"/bayestraits_N_introns.txt"),append=FALSE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(ds1[,c("code","F_intron_genes_bin","patho1")],file=paste0("./subset_",i,"/bayestraits_F_intron_genes.txt"),append=FALSE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(ds1[,c("code","LEN_introns_bin","patho1")],file=paste0("./subset_",i,"/bayestraits_LEN_introns.txt"),append=FALSE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(ds1[,c("code","LEN_exons_bin","patho1")],file=paste0("./subset_",i,"/bayestraits_LEN_exons.txt"),append=FALSE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(ds1[,c("code","LEN_intergenic_bin","patho1")],file=paste0("./subset_",i,"/bayestraits_LEN_intergenic.txt"),append=FALSE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(ds1[,c("code","tRNA_bin","patho1")],file=paste0("./subset_",i,"/bayestraits_tRNA.txt"),append=FALSE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(ds1[,c("code","pseudo_tRNA_bin","patho1")],file=paste0("./subset_",i,"/bayestraits_pseudo_tRNA.txt"),append=FALSE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  write.table(ds1[,c("code","effectors_bin","patho1")],file=paste0("./subset_",i,"/bayestraits_effectors.txt"),append=FALSE,quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
  
  
}


          