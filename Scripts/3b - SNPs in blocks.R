##### making an index of haplotype blocks
library(data.table)
library(tidyverse)
library(dplyr)
library(RColorBrewer)

colors<-brewer.pal("Paired",n=8)[c(2,6)]
multcomp<-as.numeric(read.table("Scripts/### multcomp correction value ###")[,1])

### name haplotype blocks and list snps per haplotype block
blocks<-fread("Software/XRQv1_412_239_filtered.blocks.det")
blocks$Chr_num<- as.integer(gsub("Ha412HOChr","",blocks$CHR))
blocks<- blocks %>% group_by(Chr_num) %>% mutate(hapID = paste(Chr_num,c(1:length(Chr_num)),sep="_"))
snps<-strsplit(blocks$SNPS,split="|",fixed=T)
big.list<-unlist(snps)
big.list<-data.table(SNP=big.list)
big.list$hapID<-c(rep(blocks$hapID, blocks$NSNPS))
rm(snps)

### qd solution singletons
all.snps<-fread("Software/XRQv1_412_239_filtered.map", header=F)
names(all.snps)[1:4]<-c("chr","rs","V3","ps")
all.snps$V3<-NULL

missing.snps<-all.snps[!all.snps$rs%in%big.list$SNP,]
missing.snps$Chr_num<- as.integer(gsub("Ha412HOChr","",missing.snps$chr))
missing.snps<- missing.snps %>% group_by(Chr_num) %>% mutate(hapID=paste(Chr_num,"_single",match(rs,unique(rs)),sep=""))
missing.snps<-missing.snps[,c(2,5)]
names(missing.snps)<-c("SNP","hapID")

big.list<-rbind(big.list,missing.snps)

#### read a trait to find it's haplotype blocks
thresh<-0.05/(multcomp)
suggthresh<-0.001

envs<-as.character(read.table("environments_to_run.txt")[,1])
traits<-as.character(read.table("traits_to_run.txt")[,1])

sig.blocks<-NULL #empty object to merge against
sug.blocks<-NULL
sig.snips.save<-NULL

for (i in 1:length(traits)){
  
  for(q in 1:length(envs)) {

    print(paste(traits[i],envs[q]))
    
    snips<-fread(paste("Tables/Assoc_files/", paste(traits[i],envs[q],sep="_"), ".assoc.txt", sep=""), header=T)
    

    sig.bins<-NULL
    sug.bins<-NULL

    
    tempcutoff <- as.data.frame(quantile(snips$p_wald, as.numeric(as.character(suggthresh))))[1, 1]
    sug.snips<-snips[which(snips$p_wald<tempcutoff&snips$p_wald>thresh), ]
    
    if (range(snips$p_wald)[1]<thresh){
      sig.snips<-snips[which(snips$p_wald<thresh), ]
      
      sig.snips.save<-rbind(sig.snips.save,sig.snips) ### save list of significant snps
      
      sig.bins<-merge(sig.snips,big.list,by.x="rs",by.y="SNP")
      
      sig.bins$PVE<-(2*(sig.bins$beta^2)*sig.bins$af*(1-sig.bins$af))/((2*(sig.bins$beta^2)*sig.bins$af*(1-sig.bins$af))+(2*(sig.bins$se^2)*(2*239)*sig.bins$af*(1-sig.bins$af)))
      
      sig.bins<-sig.bins %>% group_by(hapID) %>% summarise(NSNP=length(ps),beta=max(beta),min_p=min(p_wald), PVE=max(PVE))
      
      sig.bins$trait<-traits[i]
      sig.bins$env<-envs[q]
      sig.bins$pvalue<-"significant"
      print(paste(dim(sig.bins)[1],"significant haplotype blocks"))
    }
    
    
    sug.bins<-merge(sug.snips,big.list,by.x="rs",by.y="SNP")
    sug.bins$PVE<-(2*(sug.bins$beta^2)*sug.bins$af*(1-sug.bins$af))/((2*(sug.bins$beta^2)*sug.bins$af*(1-sug.bins$af))+(2*(sug.bins$se^2)*(2*239)*sug.bins$af*(1-sug.bins$af)))
    
    sug.bins<-sug.bins %>% group_by(hapID) %>% summarise(NSNP=length(ps),beta=max(beta),min_p=min(p_wald), PVE=max(PVE))
    
    if(length(sug.bins$hapID)>0){
      sug.bins$trait<-traits[i]
      sug.bins$env<-envs[q]
      sug.bins$pvalue<-"suggestive"
    }
    if(length(sug.bins$hapID)==0){
      sug.bins$trait<-NULL
      sug.bins$env<-NULL
      sug.bins$pvalue<-NULL
    }
    
    
    sug.bins<-sug.bins[!sug.bins$hapID%in%sig.bins$hapID, ]
    print(paste(dim(sug.bins)[1],"suggestive haplotype blocks"))
    
    sig.blocks<-rbind(sig.blocks,sig.bins)
    sug.blocks<-rbind(sug.blocks,sug.bins)
    
   
  }
}

sig.snips<-unique(sig.snips.save,by="rs")[,1:3]
