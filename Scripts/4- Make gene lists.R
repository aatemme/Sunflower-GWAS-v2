##############################################
#
# Author: Rishi R. Masalia, rishimasalia@gmail.com, UGA, Burke Lab, April 2017
#
# Purpose: this script is a direct follow up to SALSA by Masalia. 
# It is designed to take a "Pre_Identified_Cluster.csv" list and create a new list of genes underlying those clusters
#
# UPDATE: Stripped perl requirement (runs just in R) and started proces of simplification by Andries Temme
##########################################################################
rm(list=ls())
library(tidyverse)

genes <- 1 #buffer variable 

MoI <- read.csv(file="Tables/Bins/All-regions.csv",header=T)
MoI$X<-NULL


MoI$Run_Bin = 1:nrow(MoI)


gff3 = read.csv(file="Software/XRQ.gff3",header=F, sep="\t")
colnames(gff3) = c("Chr","Source","Feature","Start","Stop","Score","Strand","Frame","Attribute")

colnames(MoI) = c("Cat","Chr","Start","Stop", "Bin", "Run_Bin")
MoI$Start = gsub("Ha.*:","",MoI$Start,perl=T)
MoI$Stop = gsub("Ha.*:","",MoI$Stop,perl=T)

annot = read.csv(file="Software/XRQ_June2018_Annotated.csv",header=F, sep="\t")
colnames(annot)[1] = "gene_ID"
annot = annot[,c(1,3)]

biggenedf = NULL

for (m in 1:nrow(MoI)){
  if (MoI[m,2]<10){Chrom = paste("HanXRQChr0",MoI[m,2],sep="")
  } else {Chrom = paste("HanXRQChr",MoI[m,2],sep="")}
  print(Chrom)
  run_bin = MoI[m,5]
  subgff3 = gff3[which(gff3$Chr==Chrom),]
  markstart = as.numeric(as.character(MoI[m,3]))
  markstop= as.numeric(as.character(MoI[m,4]))
  
  genedf = NULL
  same = NULL
  if (markstart == markstop){same="yes"}else{same="no"}
  
  if (same == "yes"){ #same start and stop (single marker)
    ingene = NULL  
    begin = subgff3[which(subgff3$Start <= markstart),]
    arrest = subgff3[which(subgff3$Stop <= markstart),]
    nobl = as.numeric(as.character(nrow(begin)))
    noal = as.numeric(as.character(nrow(arrest)))
    
    if(nobl>noal){ingene ="yes"} else {ingene ="no"}
    
    if (ingene == "yes"){ #single marker in a gene?
      rangestart = nobl-genes
      rangestop = nobl + genes
      if (rangestart<0) { rangestart = 1}
      if (rangestop>nrow(subgff3)){rangestop = nrow(subgff3)}
      for(g in rangestart:rangestop){
        gene_ID = subgff3[g,9]
        gene_ID = gsub(";.*","",gene_ID,perl=T)
        gene_ID = gsub("ID=mRNA:","",gene_ID,perl=T)
        gene_start = subgff3[g,4]
        gene_stop = subgff3[g,5]
        
        newrow = cbind(MoI[m,c(1:4)],gene_ID, gene_start, gene_stop, run_bin)
        genedf = rbind(genedf,newrow)
      }
    } else { #single marker not in a gene (between genes)
      if (genes==1){rangestart = nobl} else { rangestart = nobl-genes }
      rangestop = nobl+genes
      if (rangestart<0) { rangestart = 1}
      if (rangestop>nrow(subgff3)){rangestop = nrow(subgff3)}
      for(g in rangestart:rangestop){
        gene_ID = subgff3[g,9]
        gene_ID = gsub(";.*","",gene_ID,perl=T)
        gene_ID = gsub("ID=mRNA:","",gene_ID,perl=T)
        gene_start = subgff3[g,4]
        gene_stop = subgff3[g,5]
        
        newrow = cbind(MoI[m,c(1:4)],gene_ID, gene_start, gene_stop, run_bin)
        genedf = rbind(genedf,newrow)
      }
    }
    
  } else { #Not same start/stop
    ingene1 = NULL  
    begin1 = subgff3[which(subgff3$Start <= markstart),]
    arrest1 = subgff3[which(subgff3$Stop <= markstart),]
    nobl1 = as.numeric(as.character(nrow(begin1)))
    noal1= as.numeric(as.character(nrow(arrest1)))
    if(nobl1>noal1){ingene1 ="yes"} else {ingene1 ="no"}
    
    ingene2 = NULL  
    begin2 = subgff3[which(subgff3$Start <= markstop),]
    arrest2 = subgff3[which(subgff3$Stop <= markstop),]
    nobl2 = as.numeric(as.character(nrow(begin2)))
    noal2= as.numeric(as.character(nrow(arrest2)))
    if(nobl2>noal2){ingene2 ="yes"} else {ingene2 ="no"}
    
    if (ingene1=="yes" && ingene2=="yes"){
      rangestart = nobl1 - genes 
      rangestop = nobl2 + genes
    }
    if (ingene1=="no" && ingene2=="yes"){
      if (genes==1){rangestart = nobl1} else { rangestart = nobl1-genes }
      rangestop = nobl2 + genes
    }
    if (ingene1=="yes" && ingene2=="no"){
      rangestart = nobl1 - genes 
      rangestop = nobl2 + genes
    }
    if (ingene1=="no" && ingene2=="no"){
      if (genes==1){rangestart = nobl1} else { rangestart = nobl1-genes }
      rangestop = nobl2 + genes
    }
    if (rangestart<0) { rangestart = 1}
    if (rangestop>nrow(subgff3)){rangestop = nrow(subgff3)}
    for(g in rangestart:rangestop){
      gene_ID = subgff3[g,9]
      gene_ID = gsub(";.*","",gene_ID,perl=T)
      gene_ID = gsub("ID=mRNA:","",gene_ID,perl=T)
      gene_start = subgff3[g,4]
      gene_stop = subgff3[g,5]
      
      newrow = cbind(MoI[m,c(1:4)],gene_ID, gene_start, gene_stop, run_bin)
      genedf = rbind(genedf,newrow)
    }
    
  } # End of Not same start/stop
  
  thing = merge(genedf,annot,by="gene_ID")
  genedf = thing[,c(2,3,4,5,1,6,7,8,9)]
  
  biggenedf = rbind(biggenedf,genedf)
} 

write.table(biggenedf,file="Tables/Genes/ListofGenes.txt",quote = F,row.names = F,col.names = T,sep="\t")

############## No Genes per Bin ##################
gene.count <- biggenedf %>% group_by(run_bin, Chr) %>% select(run_bin, gene_ID, Chr) %>%  count(run_bin)

write.table(gene.count,file="Tables/Genes/NoGenesPerBin.txt",quote = F,row.names = F,col.names = T,sep="\t")

