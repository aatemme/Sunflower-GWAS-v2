library(RColorBrewer)
library(ggpubr)
library(tidyverse)
source("Scripts/traits to haplotype blocks.R")
source("Scripts/3b- Combine snips over traits and envs.R")

### rerun haplotype analyses on just significant snps
exclude<-all.snps[!all.snps$rs%in%sig.snips$rs,]
write.table<-write.table(exclude$rs, "Tables/Blocks/snps_NOT_in_sig_blocks.txt", sep="\t", row.names=F, col.names=T, quote=F)

system("./Software/plink --tped Software/XRQv1_412_239_filtered.tped --tfam Software/XRQv1_412_239_filtered.tfam --exclude Tables/Blocks/snps_NOT_in_sig_blocks.txt --blocks 'no-pheno-req' 'no-small-max-span' --blocks-max-kb 2000000 --blocks-strong-lowci 0.7005 --out Tables/Blocks/re_sig_blocks --allow-extra-chr")


##### generate block id for snips
new.sig.blocks<-fread("Tables/Blocks/re_sig_blocks.blocks.det")

new.sig.blocks$Chr_num<- as.integer(gsub("Ha412HOChr","",new.sig.blocks$CHR))
new.sig.blocks<- new.sig.blocks %>% group_by(Chr_num) %>% mutate(hapID = paste(Chr_num,c(1:length(Chr_num)),sep="_"))
snps<-strsplit(new.sig.blocks$SNPS,split="|",fixed=T)
sig.list<-unlist(snps)
sig.list<-data.table(SNP=sig.list)
sig.list$hapID<-c(rep(new.sig.blocks$hapID, new.sig.blocks$NSNPS))
rm(snps)

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


#### calculate LD (D prime) for significant snps
system("./Software/plink --tped Software/XRQv1_412_239_filtered.tped --tfam Software/XRQv1_412_239_filtered.tfam --exclude Tables/Blocks/snps_NOT_in_sig_blocks.txt --r2 dprime yes-really --ld-window-kb 2000000 --ld-window-r2 0.0 --ld-window 1000 --out Tables/Blocks/ldtable --allow-extra-chr")

ld.table<-fread("Tables/Blocks/ldtable.ld")


### plot new blocks and ld

i<-1

chrom<-ld.table[ld.table$CHR_A==paste("Ha412HOChr",formatC(i,width=2,flag="0"), sep=""),]

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
chrom.blocks$SNP<-paste("Ha412HOChr",formatC(i,width=2,flag="0"),":",chrom.blocks$bp, sep="")



xs<-xbase+x.adjust
ys<-ybase+y.adjust
group<-c(rep(c("one"), each=4))

data<-data.frame(xs,ys,group)

a<-0.2*sin(45*pi/180)
b<-0.2*cos(45*pi/180)
c<-0.5*sin(45*pi/180)
d<-0.5*cos(45*pi/180)

xbase<-c(0.5,0.5,1.5,1.5)
x.adjust<-c(-b,-b-d,-b-d,-b)

ybase<-c(0.5,0.5,1.5,1.5)
y.adjust<-c(a,a+c,a+c,a)

plot+geom_polygon(data=data,aes(x=xs,y=ys,group=group),fill="blue")+coord_fixed()


### plot blocks

glist<-list()

for (i in 1:17) {

chrom.list<-sig.list[sig.list$chr==i, ]
chrom.list<-chrom.list %>% gather(key="hap_type", value="block_id", hapID, sigblock_hapID)
chrom.list$type_id<-paste(chrom.list$hap_type,chrom.list$block_id,sep="-")


colours<-rep(c(brewer.pal(8,"Dark2"))[-7],length(unique(chrom.list$type_id)))

plot<-ggplot(data= chrom.list, aes(x=fct_reorder(SNP,ps),y=hap_type,fill=fct_inorder(type_id)))+
  geom_tile()+
  scale_fill_manual(values=colours)+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90, vjust = 0.5,hjust=1),axis.title.x = element_blank())+
  ggtitle(paste("chrom",i))

glist[[i]]<-plot

}

multi.page <- ggpubr::ggarrange(plotlist = glist, nrow = 3, ncol = 1)
ggpubr::ggexport(multi.page, filename = "Plots/Colocalization/oldblocks vs newblocks.pdf")
