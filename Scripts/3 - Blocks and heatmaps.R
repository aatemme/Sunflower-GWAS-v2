library(RColorBrewer)
library(ggpubr)
library(tidyverse)
library(grid)
library(ggrepel)
source("Scripts/3b - SNPs in blocks.R")

#### read in preferences
prefs<-read.table("Scripts/### Preferences ###",header=F,sep="=",skip=1)
  SNPset<-as.character(prefs[2,2])
  pheno.name<-as.character(prefs[1,2])
  multcomp<-as.numeric(as.character(prefs[3,2]))

### rerun haplotype analyses on just significant snps

exclude<-big.list[!big.list$SNP%in%sig.snips$rs,]

write.table<-write.table(exclude$SNP, "Tables/Blocks/snps_NOT_in_sig_blocks.txt", sep="\t", row.names=F, col.names=T, quote=F)

system(paste("./Software/plink --tped Software/",SNPset,".tped --tfam Software/",SNPset,".tfam --exclude Tables/Blocks/snps_NOT_in_sig_blocks.txt --blocks 'no-pheno-req' 'no-small-max-span' --blocks-max-kb 2000000 --blocks-strong-lowci 0.7005 --out Tables/Blocks/re_sig_blocks --allow-extra-chr --blocks-inform-frac 0.9",sep=""))


##### generate block id for snips
new.sig.blocks<-fread("Tables/Blocks/re_sig_blocks.blocks.det")

if (dim(new.sig.blocks)[1]>0) {
new.sig.blocks$Chr_num<- as.integer(gsub("Ha412HOChr","",new.sig.blocks$CHR))
new.sig.blocks<- new.sig.blocks %>% group_by(Chr_num) %>% mutate(hapID = paste(Chr_num,c(1:length(Chr_num)),sep="_"))
snps<-strsplit(new.sig.blocks$SNPS,split="|",fixed=T)
sig.list<-unlist(snps)
sig.list<-data.table(SNP=sig.list)
sig.list$hapID<-c(rep(new.sig.blocks$hapID, new.sig.blocks$NSNPS))
rm(snps)
} else {sig.list<-data.frame(SNP=c(),hapID=c())} ## for when there are zero blocks

##### add in in singletons for significant snps
missing.snps<-sig.snips[!sig.snips$rs%in%sig.list$SNP, c(1:3) ]
missing.snps$Chr_num<- as.integer(gsub("Ha412HOChr","",missing.snps$chr))
missing.snps<- missing.snps %>% group_by(Chr_num) %>% mutate(hapID=paste(Chr_num,"_single-sig",match(rs,unique(rs)),sep=""))
missing.snps<-missing.snps[,c(2,5)]
names(missing.snps)<-c("SNP","hapID")

sig.list<-rbind(sig.list,missing.snps)


### split SNP name into info
sig.list$chr <- gsub("Ha412HOChr","",sig.list$SNP)
sig.list$chr <- as.integer(gsub(":.*","",sig.list$chr))
sig.list$ps <- as.integer(gsub(".*:","",sig.list$SNP))

sig.list<-sig.list[order(chr,ps),]


#### add big haplotype block id's
sig.list$sigblock_hapID<-sig.list$hapID
sig.list$hapID<-NULL

sig.list$hapID<-big.list$hapID[match(sig.list$SNP,big.list$SNP)]

sig.list<-sig.list[order(sig.list$chr,sig.list$ps), ] ## order

sig.list$sigblock_hapID<-fct_inorder(sig.list$sigblock_hapID)


#### flag haplotype blocks that have been split in the new blocking instead of combined fully
split.blocks<-sig.list %>% group_by(hapID) %>% summarise(blocks=length(unique(sigblock_hapID))) %>% filter(blocks>1)

      ### keep full genome haplotype block when the block isn't fully merged with another
      sig.list$sigblock_hapID<-as.character(sig.list$sigblock_hapID)
      sig.list$sigblock_hapID[sig.list$hapID%in%split.blocks$hapID]<-paste(as.character(sig.list$hapID)[sig.list$hapID%in%split.blocks$hapID],"kept")

sig.list$sigblock_hapID<-fct_inorder(sig.list$sigblock_hapID)

#### pretty new block name
sig.list<-sig.list %>% group_by(chr) %>%
  mutate(region=paste(formatC(as.numeric(chr),width=2, flag="0"), 
                      formatC(as.numeric(factor(rank(match(sigblock_hapID,levels(sigblock_hapID))))),width=2, flag="0"),sep="-")) # super janky ordering

#### add new (combined) block name to significant blocks data
sighap_to_genomehap<-data.frame(genome.hap=unique(sig.list$hapID))
sighap_to_genomehap$sig.hap<-sig.list$sigblock_hapID[match(sighap_to_genomehap$genome.hap,sig.list$hapID)]
sighap_to_genomehap$colocate.region<-sig.list$region[match(sighap_to_genomehap$genome.hap,sig.list$hapID)]

### save some of the blocks objects for later
write.table<-write.table(sig.blocks, "Tables/Blocks/traits_to_genomeblocks_signif.txt", sep="\t", row.names=F, col.names=T)
write.table<-write.table(sug.blocks, "Tables/Blocks/traits_to_genomeblocks_sugest.txt", sep="\t", row.names=F, col.names=T)
write.table<-write.table(sig.list, "Tables/Blocks/sigsnips_to_genomeblocks.txt", sep="\t", row.names=F, col.names=T)
write.table<-write.table(sighap_to_genomehap, "Tables/Blocks/condensed_genome_blocks.txt", sep="\t", row.names=F, col.names=T)

#### calculate LD (D prime) for significant snps

system(paste("./Software/plink --tped Software/",SNPset,".tped --tfam Software/",SNPset,".tfam --exclude Tables/Blocks/snps_NOT_in_sig_blocks.txt --r2 dprime yes-really --ld-window-kb 2000000 --ld-window-r2 0.0 --ld-window 1000 --out Tables/Blocks/ldtable --allow-extra-chr",sep=""))

ld.table<-fread("Tables/Blocks/ldtable.ld")


### plot new blocks and ld

for (i in 1:length(unique(ld.table$CHR_A))) {
chrom<-ld.table[ld.table$CHR_A==unique(ld.table$CHR_A)[i],]

Chr.num<-as.numeric(gsub("Ha412HOChr","",unique(ld.table$CHR_A)[i]))

chrom.snps<-unique(c(as.character(chrom$BP_A),as.character(chrom$BP_B)))
nsnps<-length(chrom.snps)[1]

chrom.snps<-data.table(BP_A=as.numeric(chrom.snps),BP_B=as.numeric(chrom.snps))

chrom<-rbind(chrom,chrom.snps,fill=T)

chrom$bp_A<-fct_reorder(factor(chrom$BP_A),chrom$BP_A)
chrom$bp_B<-fct_reorder(factor(chrom$BP_B),chrom$BP_B)


## draw LD plot
plot<-ggplot(data=chrom,aes(y=bp_A, x=bp_B,fill=R2))+geom_tile()+scale_fill_gradient(low = "white", high = "#E64A19",na.value = "black")+
  annotate(geom="polygon",x=c(0,0,nsnps+1), y=c(0,nsnps+1,nsnps+1),fill="white")+scale_y_discrete(position = "right")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1),axis.text.y = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),axis.title.y = element_blank())


## get blocks info for the snps
chrom.blocks<-data.frame(bp=unique(chrom$BP_A))
chrom.blocks$SNP<-paste("Ha412HOChr",formatC(Chr.num,width=2,flag="0"),":",chrom.blocks$bp, sep="")
chrom.blocks$big.hap<-fct_inorder(sig.list$hapID[match(chrom.blocks$SNP,sig.list$SNP)])
chrom.blocks$sig.hap<-fct_inorder(as.character(sig.list$sigblock_hapID[match(chrom.blocks$SNP,sig.list$SNP)]))
chrom.blocks$colocate.region<-sig.list$region[match(chrom.blocks$SNP,sig.list$SNP)]

#### add block info to the ld plot
colours<-rep(c(brewer.pal(8,"Dark2"))[-7],100)
chrom.blocks$big.hap.colors<-colours[as.numeric(chrom.blocks$big.hap)]
chrom.blocks$sig.hap.colors<-colours[as.numeric(chrom.blocks$sig.hap)]

a<-0.2*(length(chrom.blocks$SNP)/10)*sin(45*pi/180)
b<-0.2*(length(chrom.blocks$SNP)/10)*cos(45*pi/180)
c<-0.5*(length(chrom.blocks$SNP)/10)*sin(45*pi/180)
d<-0.5*(length(chrom.blocks$SNP)/10)*cos(45*pi/180)

xbase<-c(rep(c(0:(length(chrom.blocks$SNP)-1)),each=4))+c(0.5,0.5,1.5,1.5)
x.adjust<-c(-b,-b-d,-b-d,-b)
x2.adjust<-c(-b-d,-b-d-d,-b-d-d,-b-d)

ybase<-c(rep(c(0:(length(chrom.blocks$SNP)-1)),each=4))+c(0.5,0.5,1.5,1.5)
y.adjust<-c(a,a+c,a+c,a)
y2.adjust<-c(a+c,a+c+c,a+c+c,a+c)

xs<-xbase+x.adjust
xs2<-xbase+x2.adjust

ys<-ybase+y.adjust
ys2<-ybase+y2.adjust
group<-c(rep(c(0:(length(chrom.blocks$SNP)-1)),each=4))

big.hap<-data.frame(xs,ys,group)
sig.hap<-data.frame(xs2,ys2,group)

sig.hap$colocate.region<-rep(chrom.blocks$colocate.region,each=4)

test<-sig.hap %>% group_by(colocate.region) %>% summarize (x=mean(xs2),y=mean(ys2))

hap.plot<-plot+geom_polygon(data=big.hap,aes(x=xs,y=ys,group=group),fill=rep(chrom.blocks$big.hap.colors,each=4))+
      geom_polygon(data=sig.hap,aes(x=xs2,y=ys2,group=group),fill=rep(chrom.blocks$sig.hap.colors,each=4))+
  geom_segment(x=xs2[1],y=ys2[1],xend=xs2[length(xs2)],yend=ys2[length(ys2)],col="black")+
coord_fixed(xlim=c(0.5,(length(chrom.blocks$SNP)+0.5)*1.05),ylim=c(0.5-(length(chrom.blocks$SNP)+0.5)*0.05,length(chrom.blocks$SNP)+0.5),clip="off",expand=0)+
  theme(legend.position="none")+
  annotate("text", x=xs[1],y=ys[1],label="Genome ", angle=45,hjust=1,vjust=0)+
  annotate("text", x=xs2[1],y=ys2[1],label="Significant ", angle=45,hjust=1,vjust=0)
  # annotate("text",x=test$x,y=test$y,label=test$colocate.region, angle=45,hjust=0.5,vjust=0.5)


hap.plot<-hap.plot+geom_text_repel(data=test,aes(x=x,y=y,label=colocate.region,group=NULL,fill=NULL),
                         nudge_y=8*a,
                         nudge_x=-b,
                         direction="y",
                         angle=45)

legend<-get_legend(plot)

pdf(paste("Plots/Colocalization/Chromosome-",Chr.num,".pdf",sep=""),height=7.5,width=10.5)

# png(paste("Plots/Colocalization/Chromosome-",i,".png",sep=""),height=750,width=1050)

grid.newpage()

pushViewport(viewport(name = "rotate", angle = -45, clip = "off", width = 0.9, height = 0.9))
print(hap.plot, vp = "rotate")

vp<-viewport(x=0.28,y=0.8,width=0.4,height = 0.1)
pushViewport(vp)
grid.text(as.character(paste("chromosome:",Chr.num)), 0.2, 0.2,gp=gpar(cex=1.5))
popViewport(1)

vp<-viewport(x=0.6,y=0.85,width=0.1,height = 0.1)
pushViewport(vp)
grid.draw(legend)
popViewport(1)
# popViewport(1)
dev.off()

}

