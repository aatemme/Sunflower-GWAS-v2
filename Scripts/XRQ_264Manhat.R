#########################################################
### Burrito {Manhattanizer}
#########################################################
#
# Author: Sariel Hubner, sariel.hubner@botany.ubc.ca
# Ammended by: Rishi R. Masalia, rishimasalia@gmail.com, Jan 2016 - April 2017, Burke Lab, UGA
# Please don't distribute
#
# Summary: Burrito is a wrapper for the Manhattanizer (augmentation of Sariel's original Helimap pipeline)
# It treats traits and environments individually, and produces a Manhattan plot & list of significant SNPs.
# EMMAX is run with kinship and structure (PCA) to get associations
#
# 
#
#########################################################
### Installing and loading required packages
#########################################################

if (!require("qqman")) {
  install.packages("qqman", dependencies = TRUE)
  library(qqman)
}
if (!require("CMplot")) {
  install.packages("CMplot", dependencies = TRUE)
  library(CMplot)
}

if (!require("SNPRelate")) {
  source("https://bioconductor.org/biocLite.R")
  biocLite("SNPRelate")
  library(SNPRelate)
}

if (!require("lme4")) {
  install.packages("lme4", dependencies = TRUE)
  library(lme4)
}

if (!require("gplots")) {
  install.packages("gplots", dependencies = TRUE)
  library(gplots)
}

if (!require("metafor")) {
  install.packages("metafor", dependencies = TRUE)
  library(metafor)
}

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer", dependencies = TRUE)
  library(RColorBrewer)
}

library(tools)


#### getting this script working without burrit wrapper

envs<-as.character(read.table("Environments.list")[,1])
traits<-as.character(read.table("Traits.list")[,1])



#########################################################
### Read and prepare files
#########################################################

# setwd("/home/atemme/Desktop/Seedmorph/EnchiladaSuite/Burrito") ##### change this to match your path to Enchilada

##########################################################33
# config = read.table(file="Asada.config")

trait = traits[1]
envNames <- envs[1]
envFiles <- paste(trait,"_",envNames,".txt",sep="")
num_env <-1
num_rep <-1
prefix <- "XRQ_264"


#########read phenotypes and split to have data frame for each env and rep:

op <- options(stringsAsFactors=FALSE)
sam<-read.table("Software/XRQ_264.tfam")
SAM<-sam[,1]
options(op)

#########################################################
### Association
#########################################################
# num_env<-as.numeric(as.character(num_env))
# num_rep<-as.numeric(as.character(num_rep))
# print(num_env)
# print(num_rep)



Pheno_in<-read.table(paste("Phenotype data/",envFiles[1],sep=""),header=TRUE)

Pheno_in[Pheno_in=="Na"]<-"NA"

PhenFile<-merge(SAM, Pheno_in, by.x=1,by.y=1,all.x=TRUE)


	options(op)
	write.table(data.frame(PhenFile[,1],PhenFile[,1], PhenFile[,2]),paste("Software/Temporary files/",trait,"_",envNames[1],".pheno",sep=""),row.names = FALSE,col.names = FALSE, quote = FALSE)
	# if (model == "K" || model == "k"){;
	  print ("P+K")
	  EMMAX.pk<-paste("./Software/emmax-intel64 -v -d 10 -t ",paste("./Software/",prefix,sep=""), " -p ", paste("./Software/",trait,"_",envNames[1],".pheno",sep=""), " -c ", paste("./Software/",prefix,".PCA_EV",sep=""), " -k ", paste("./Software/",prefix,".aIBS.kinf",sep=""), " -o ", paste(trait,"_",envNames[1],sep=""), sep="")
	  system (EMMAX.pk)

	system (EMMAX.env1)

  move<-paste("mv ./",paste(trait,"_",envNames[1],".ps",sep=""),paste(" ./Tables/PSfiles/",sep=""),sep="")

  system(move)


system(paste("rm ",paste(trait,"_",envNames[1],sep=""),".log",sep=""))
system(paste("rm ",paste(trait,"_",envNames[1],sep=""),".reml",sep=""))
system(paste("rm Software/",paste(trait,"_",envNames[1],sep=""),".pheno",sep=""))





