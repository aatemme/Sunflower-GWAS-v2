##### making an index of haplotype blocks
library(data.table)
library(tidyverse)
library(dplyr)
library(RColorBrewer)

colors<-brewer.pal("Paired",n=8)[c(2,6)]


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
all.snps<-fread("Tables/Assoc_files/Calcium_logdiff.assoc.txt", header=T)

missing.snps<-all.snps[!all.snps$rs%in%big.list$SNP, c(1:3) ]
missing.snps$Chr_num<- as.integer(gsub("Ha412HOChr","",missing.snps$chr))
missing.snps<- missing.snps %>% group_by(Chr_num) %>% mutate(hapID=paste(Chr_num,"_single",match(rs,unique(rs)),sep=""))
missing.snps<-missing.snps[,c(2,5)]
names(missing.snps)<-c("SNP","hapID")

big.list<-rbind(big.list,missing.snps)

#### read a trait to find it's haplotype blocks
thresh<-0.05/(31733+15809)
suggthresh<-0.001

envs<-as.character(read.table("environments_to_run.txt")[,1])
traits<-as.character(read.table("traits_to_run.txt")[,1])

sig.blocks<-NULL #empty object to merge against
sug.blocks<-NULL

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
      
      sig.bins<-merge(sig.snips,big.list,by.x="rs",by.y="SNP")
      
      sig.bins<-sig.bins %>% group_by(hapID) %>% summarise(NSNP=length(ps),beta=median(beta),min_p=min(p_wald))
      
      sig.bins$trait<-traits[i]
      sig.bins$env<-envs[q]
      sig.bins$pvalue<-"significant"
      print(paste(dim(sig.bins)[1],"significant haplotype blocks"))
    }
    
    
    sug.bins<-merge(sug.snips,big.list,by.x="rs",by.y="SNP")
    sug.bins<-sug.bins %>% group_by(hapID) %>% summarise(NSNP=length(ps),beta=median(beta),min_p=min(p_wald))
    
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

# exclude<-big.list[!big.list$hapID%in%unique(c(sig.blocks$hapID)),]
# 
# write.table<-write.table(exclude$SNP, "Tables/snps_NOT_in_sig_blocks.txt", sep="\t", row.names=F, col.names=T, quote=F)
# 
# 
# system("./Software/plink --tped Software/XRQv1_412_239_filtered.tped --tfam Software/XRQv1_412_239_filtered.tfam --exclude Tables/snps_NOT_in_sig_blocks.txt --blocks 'no-pheno-req' 'no-small-max-span' --blocks-max-kb 40000 --blocks-strong-lowci 0.7005 --out Tables/re_sig_blocks --allow-extra-chr")
# 
# 
# new.sig.blocks<-fread("Tables/re_sig_blocks.blocks.det")


# colocate<-rbind(sig.blocks,sug.blocks[sug.blocks$hapID%in%sig.blocks$hapID,])
# 
# colocate$trait_env<-paste(colocate$trait,colocate$env,sep="_")
# 
# traits.per.block<-colocate %>% group_by(hapID) %>% summarise(trait_num=length(trait_env))
# 
# single.trait.blocks<-colocate[colocate$hapID%in%traits.per.block$hapID[traits.per.block$trait_num==1],]
# 
# colocate<-colocate[!colocate$hapID%in%traits.per.block$hapID[traits.per.block$trait_num==1],]
# 
# colocate<-colocate %>% separate(hapID, sep= "_", c("chromosome","blocknum"),remove=F) %>% 
#                         arrange(chromosome, blocknum)
# 
# colocate<- colocate %>% group_by(chromosome) %>% 
#                         mutate(blockseq=match(blocknum,unique(blocknum)),beta.sign=sign(beta)) %>%
#                         mutate(region=paste(formatC(as.numeric(chromosome),width=2, flag="0"), formatC(blockseq,width=2, flag="0"),sep="-"))
# 
# chrom.borders<-colocate %>% group_by(chromosome)%>% summarise(bin=n_distinct(blockseq))
# chrom.borders<-cumsum(chrom.borders$bin)
# chrom.borders<-chrom.borders+0.5
# chrom.borders<-chrom.borders[1:length(chrom.borders)-1]
# 
# 
# 
# 
# plot.data<-colocate[colocate$env=="water",]
# 
# chrom.borders<-plot.data %>% group_by(chromosome)%>% summarise(bin=length(unique(region))) %>% arrange(as.integer(chromosome))
# chrom.borders<-cumsum(chrom.borders$bin)
# chrom.borders<-chrom.borders+0.5
# chrom.borders<-chrom.borders[1:length(chrom.borders)-1]
# 
# baseplot<-ggplot(plot.data,aes(x=region,y=trait,fill=as.factor(beta.sign)))
# 
# baseplot+geom_vline(xintercept=c(1:length(plot.data$region)),colour="darkgrey",linetype=3)+
#   geom_vline(xintercept=chrom.borders,colour="black")+
#   geom_tile(fill="white")+
#   geom_tile(aes(alpha=pvalue),colour="black")+
#   theme_minimal()+
#   theme(axis.text.y = element_text(hjust = 0))+
#   scale_fill_manual(values=c(colors[1],colors[2]))+
#   scale_alpha_manual(values=c(1,0.1))+
#   scale_x_discrete(drop=F)+
#   theme_classic()+
#   theme(axis.title.y=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5,hjust=1))+
#   ggtitle("control")+theme(legend.position = "none")+theme(axis.title.x=element_blank())
# 
# 

