library(RColorBrewer)
source("Scripts/traits to haplotype blocks.R")
source("Scripts/3b- Combine snips over traits and envs.R")




### rerun haplotype analyses on just significant snps
exclude<-all.snps[!all.snps$rs%in%sig.snips$rs,]
write.table<-write.table(exclude$rs, "Tables/snps_NOT_in_sig_blocks.txt", sep="\t", row.names=F, col.names=T, quote=F)

system("./Software/plink --tped Software/XRQv1_412_239_filtered.tped --tfam Software/XRQv1_412_239_filtered.tfam --exclude Tables/snps_NOT_in_sig_blocks.txt --blocks 'no-pheno-req' 'no-small-max-span' --blocks-max-kb 2000000 --blocks-strong-lowci 0.7005 --out Tables/re_sig_blocks --allow-extra-chr")

##### generate block id for snips
new.sig.blocks<-fread("Tables/re_sig_blocks.blocks.det")

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

### plot blocks

chrom.list<-sig.list[sig.list$chr==10, ]
chrom.list<-chrom.list %>% gather(key="hap_type", value="block_id", hapID, sigblock_hapID)
chrom.list$type_id<-paste(chrom.list$hap_type,chrom.list$block_id,sep="-")


colours<-rep(c(brewer.pal(8,"Dark2"))[-7],length(unique(chrom.list$type_id)))

ggplot(data= chrom.list, aes(x=fct_reorder(SNP,ps),y=hap_type,fill=fct_inorder(type_id)))+
  geom_tile()+
  scale_fill_manual(values=colours)+
  theme(legend.position = "none",axis.text.x = element_text(angle = 90, vjust = 0.5,hjust=1))



