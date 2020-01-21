
#### visualizing significant blocks on the haplotype map
###### tweaking haplotype blocks estimation
library(data.table)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(RColorBrewer)
library(scales)

#### read in preferences
prefs<-read.table("Scripts/### Preferences ###",header=F,sep="=",skip=1)
    SNPset<-as.character(prefs[2,2])
    pheno.name<-as.character(prefs[1,2])
    multcomp<-as.numeric(as.character(prefs[3,2]))

####


#### name blocks/snps <- same process as other steps where blocks are named == same names
      blocks<-fread(paste("Software/",SNPset,".blocks.det",sep=""))
      
      blocks$Chr_num<- as.integer(gsub("Ha412HOChr","",blocks$CHR))
      blocks<- blocks %>% group_by(Chr_num) %>% mutate(hapID = paste(Chr_num,c(1:length(Chr_num)),sep="_"))
      snps<-strsplit(blocks$SNPS,split="|",fixed=T)
      big.list<-unlist(snps)
      big.list<-data.table(SNP=big.list)
      big.list$hapID<-c(rep(blocks$hapID, blocks$NSNPS))
      rm(snps)
      
      ### qd solution singletons
      all.snps<-fread(paste("Software/",SNPset,".map",sep=""), header=F)
      names(all.snps)[1:4]<-c("chr","rs","V3","ps")
      all.snps$V3<-NULL
      
      missing.snps<-all.snps[!all.snps$rs%in%big.list$SNP,]
      missing.snps$Chr_num<- as.integer(gsub("Ha412HOChr","",missing.snps$chr))
      missing.snps<- missing.snps %>% group_by(Chr_num) %>% mutate(hapID=paste(Chr_num,"_single",match(rs,unique(rs)),sep=""))
      missing.snps<-missing.snps[,c(2,5)]
      names(missing.snps)<-c("SNP","hapID")
      missing.snps<-separate(data=missing.snps,col=SNP, c("CHR", "BP1"),sep=":",remove=F)
      missing.snps$BP1<-as.numeric(missing.snps$BP1)
      missing.snps$BP2<-missing.snps$BP1+1
      
##### Identify blocks with significant snps
      
      #######
      sig.list<-read.table("Tables/Blocks/sigsnips_to_genomeblocks.txt",header=T)
      genemap<-read.table("Tables/Blocks/condensed_genome_blocks.txt",header=T)
      # colocate<-read.table("Tables/Blocks/colocate_table.txt")
      genemap$colocate.block<-genemap$colocate.region #rename
      
      genemap$start<-blocks$BP1[match(genemap$genome.hap,blocks$hapID)]
      genemap$stop<-blocks$BP2[match(genemap$genome.hap,blocks$hapID)]
      genemap$chr<-blocks$CHR[match(genemap$genome.hap,blocks$hapID)]
      genemap$block.type<-"region"
      
      single.snps<-sig.list$SNP[match(genemap$genome.hap[is.na(genemap$start)],sig.list$hapID)]
      single.snps.ps<-sig.list$ps[match(genemap$genome.hap[is.na(genemap$start)],sig.list$hapID)]
      single.snps.chr<-paste("Ha412HOChr",formatC(sig.list$chr[match(genemap$genome.hap[is.na(genemap$start)],sig.list$hapID)],width =2,flag="0"),sep="")
      
      genemap$block.type[is.na(genemap$start)]<-"single"
      genemap$start[is.na(genemap$start)]<-single.snps.ps
      genemap$stop[is.na(genemap$stop)]<-single.snps.ps
      genemap$chr[is.na(genemap$chr)]<-single.snps.chr
      genemap$pointpos<-(genemap$start+genemap$stop)/2
      genemap$CHR<-genemap$chr
      

#### block sizes point 9
nrcol<-4
blocks$blockcol<-rep(1:nrcol,length=nrow(blocks))

plotbase<-ggplot(data=blocks,aes(x=BP1,y=1))
blockmap<-plotbase+theme_minimal()+
  geom_rect(data=blocks,xmin=blocks$BP1,xmax=blocks$BP2,ymin=0,ymax=1,aes(fill=as.factor(blockcol)))+
  scale_fill_manual(values=brewer.pal(nrcol,"Dark2"))+
  # geom_segment(x=blocks$BP1,xend=blocks$BP2,y=0,yend=0,col="black")+
  geom_segment(data=missing.snps,x=missing.snps$BP1, xend=missing.snps$BP2+10000,y=0,yend=0,col="black")+
  # geom_point(data=missing.snps,x=missing.snps$BP1,y=1,col="grey50",size=0.2)+
  geom_point(data=genemap, x=genemap$pointpos,y=0,col="gray50",size=0.5)+coord_cartesian(ylim=c(0,1))+
  facet_wrap(~CHR,scales="free_y",ncol=3,nrow=6,strip.position = "left")+theme(legend.position="none")+
  scale_x_continuous(labels = comma,breaks=c(0,50e6,100e6,150e6,200e6))+
  theme(axis.text.x = element_text(angle=20,hjust=1),
        axis.text.y = element_blank(),axis.title.y = element_blank())+xlab("Base position")

ggsave("Plots/Colocalization/regions_on_blockmap.pdf",blockmap, width=14.24, height=7.6)
