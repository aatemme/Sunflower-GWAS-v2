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
# color1 = as.character(config[7,1])
# color2 = as.character(config[8,1])
model <- "K" #can also run with fastSTRUCTURE groups
comparison_name = envNames[1] 

#########read phenotypes and split to have data frame for each env and rep:

op <- options(stringsAsFactors=FALSE)
sam<-read.table("Software/Genomics data/XRQ_264.tfam")
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
	  EMMAX.pk<-paste("./Software/emmax-intel64 -v -d 10 -t ",prefix, " -p ", paste("./Software/",trait,"_",envNames[1],".pheno",sep=""), " -c ", paste("./Software/",prefix,".PCA_EV",sep=""), " -k ", paste("./Software/",prefix,".aIBS.kinf",sep=""), " -o ", paste(trait,"_",envNames[1],sep=""), sep="")
	  system (EMMAX.pk)
	# }else {
	# print ("Q+K model being run")
	# EMMAX.env1<-paste("./bin/emmax-intel64 -v -d 10 -t ",
	#                   prefix, 
	#                   " -p ",
	#                   paste(trait,"_",envNames[i],".pheno",sep=""),
	#                   " -c FS_K5.structure",
	#                   " -k ", paste(prefix,".aIBS.kinf",sep=""),
	#                   " -o ", paste(trait,"_",envNames[i],sep=""), sep="")
	system (EMMAX.env1)

  


cat ("Reading the map.......................")
snp.map = read.table(file="XRQ_264.map",header= T) 

numLoci = 649181 
InferrMeff = 92747 
adj = (0.05/92747)
cat("Total number of SNPs is: ", numLoci, "\n")
cat("Inferred Meff is: ", InferrMeff, "\n")
#adj<-0.05/sum(simpleMeff)
adjP<-format(adj, scientific = FALSE)

cat("Corrected threshold is: ", adjP, "\n")

#########################################################
### Summarize results
#########################################################
  cat("summarizing results.............\n")
  
  reps <- Sys.glob("*.ps")
  
  out.file<-data.frame()
  for (r in reps)
  {  
    Fname<- file_path_sans_ext(r)
    ps <- read.delim(r,header=FALSE)
    Sig<-data.frame()
    inFile<-read.delim(r,header=FALSE)
    print(r)
    Sig<-inFile[inFile[,4]<=adj,] #adj
    if (nrow(Sig)>0) {
      ifelse (nrow(Sig)>0, Sig[,5]<-Fname,NA)
      out.file <- rbind(out.file, Sig)
    }
  	
    ######## If there are NO significant Values in the current looped environment######################
     if (nrow(Sig)==0) { next }  
    ##################################################################

  }
if (nrow(Sig)==0){
  Wording = "There are no significant SNPs in this trait in this environment"
  print (Wording)
  cat(Wording, file=(paste(trait,"_",comparison_name,".sigsnps",sep="")), sep="\t", fill=FALSE, append=FALSE)
  # reps <- Sys.glob("*.ps")
  # ColSet<-c(color1,color2)
  # 
  # pdf(paste(comparison_name,"_ManhattanPlot_",trait,".pdf",sep=""))
  # if(num_env>=num_rep) {
  #   par(mfrow=(c(num_env,num_rep))) } else {par(mfrow=(c(num_rep,num_env)))}
  # for(r in reps) {
  #   Fname<- file_path_sans_ext(r)
  #   ps <- read.delim(r,header=FALSE)
  #   pvalue = ps$V4
  #   ytop = -log10(min(pvalue))
  #   if (ytop > 10){ ytop = ytop+2} else { ytop = 12}
  #   ps1<- merge(snp.map,ps[,c(1,4)],by.x=1,by.y=1)
  #   colnames(ps1) = c("SNP_ID", "Chr", "Pos", "V4")
  #   print(r)
  #   NNN<-strsplit(Fname,"_")[[1]]
  #   envN<-NNN[2]
  #  #Print with qfdr:
  #   manhattan(ps1[ps1$V4<0.01,],chr = "Chr", bp = "Pos", p = "V4", snp = "SNP_ID",main=paste(trait," ",comparison_name,sep=""),col=ColSet,suggestiveline = FALSE, genomewideline=-log10(adj), cex=0.6, cex.axis=0.8,ylim=c(2,ytop))
  #   }
  #   dev.off()
  #   
  #   cat("starting cleanup.............\n")
  #   system(paste("rm *.log",sep=""))
  #   system(paste("rm *.reml",sep=""))
  #   system(paste("rm *.pheno",sep=""))
  #   system(paste("rm *.txt",sep=""))
  #   cat("All Done!.............\n")
  quit()
} ##### This is a critical break; if there are 0 significant SNPs in this trait env , so it breaks

SigSNPs<-out.file
colnames(SigSNPs)<-c("SNP_ID","beta","SE_beta","Pvalue","DataSet")

phenotypes<-Sys.glob("*.pheno")
phenFile<-data.frame()
for (phen in phenotypes)
{
  Pname<- file_path_sans_ext(phen)
  print(phen)
  Var<-data.frame(Pname,"NA")
  trFile <- read.table(phen,header=FALSE)
  Var[,2]<-var(trFile[,3],na.rm=TRUE)
 # ifelse (is.na(Var), Var[,2]<-phen,NA)
  phenFile <- rbind(phenFile, Var)
}

Var.file<-phenFile
colnames(Var.file)<-c("DataSet","Phenotypic_Variance")
SumFile<-merge(SigSNPs,Var.file,by.x=5,by.y=1,all.x=TRUE)


SumFile$Phenotypic_Variance<-as.numeric(SumFile$Phenotypic_Variance)
SumFile$SE_beta<-as.numeric(SumFile$SE_beta)
PVE<-100*(SumFile$SE_beta^2/SumFile$Phenotypic_Variance)
Summary<-data.frame(SumFile[,1:6],PVE)

colnames(snp.map)<-c("SNP_ID","Chr","Pos")
Sum2<-merge(snp.map,Summary,by.x=1,by.y=2)
Sum2$Env = 0
for(val in 1:nrow(Sum2)){
 env_name = sub("Ha412HO_April2017:","",Sum2[val,4])
 env_name = sub(":.*","",env_name)
 Sum2$Env[val] <- env_name
 #Sum2$Chr = gsub("nXRQChr","",Sum2$Chr)
}
#Sum2$Chr = as.numeric(as.character(Sum2$Chr))

write.table(Sum2, paste(trait,"_",comparison_name,".sigsnps",sep=""),row.names = FALSE, quote = FALSE,na="NA",sep="\t")

#########################################################
### Plot results
#########################################################
# 
# cat ("Running Manhattan Plots.........\n")
# reps <- Sys.glob("*.ps")
# ColSet<-c(color1,color2)
# 
# pdf(paste(comparison_name,"_ManhattanPlot_",trait,".pdf",sep=""))
# 
# for(r in reps) {
#   Fname<- file_path_sans_ext(r)
#   ps <- read.delim(r,header=FALSE)
#   pvalue = ps$V4
#   ytop = -log10(min(pvalue))
#   if (ytop > 10){ ytop = ytop+2} else { ytop = 12}
#   ps1<- merge(snp.map,ps[,c(1,4)],by.x=1,by.y=1)
#   print(r)
#   NNN<-strsplit(Fname,":")[[1]]
#   envN<-NNN[2]
#   repN<-NNN[3]
#   
#   #Print no qfdr:
#   manhattan(ps1[ps1$V4<0.01,],chr = "Chr", bp = "Pos", p = "V4", snp = "SNP_ID",main=paste(trait," ",comparison_name,sep=""),col=ColSet,suggestiveline = FALSE,genomewideline=-log10(adj), cex=0.6, cex.axis=0.8,ylim=c(2,ytop))
#   
#   #Print with qfdr:
#   #manhattan(ps1[ps1$V4<0.01,],chr = "Chr", bp = "Pos", p = "V4", snp = "SNP_ID",main=paste(trait," ",comparison_name,sep=""),col=ColSet,suggestiveline=-log10(qfdr),genomewideline=-log10(adj), cex=0.6, cex.axis=0.8,ylim=c(2,ytop))
#   
#   #Print w no qfdr:
#   #manhattan(ps1[ps1$V4<0.01,],chr = "Chr", bp = "Pos", p = "V4", snp = "SNP_ID",main=paste(trait," ",comparison_name,sep=""),col=ColSet,genomewideline=-log10(adj), cex=0.6, cex.axis=0.8,ylim=c(2,ytop))
#   
# }
# dev.off()

#########################################################
### Cleanup 
#########################################################
cat("starting cleanup.............\n")
system(paste("rm *.log",sep=""))
system(paste("rm *.reml",sep=""))
system(paste("rm *.pheno",sep=""))
system(paste("rm *.txt",sep=""))
cat("All Done!.............\n")



