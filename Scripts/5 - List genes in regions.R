### new genelist finder
library(tidyverse)
library(data.table)
library(urltools) #for parsing of gene product special characters

#### read in preferences
prefs<-read.table("Scripts/### Preferences ###",header=F,sep="=",skip=1)
  SNPset<-as.character(prefs[2,2])
  pheno.name<-as.character(prefs[1,2])
  multcomp<-as.numeric(as.character(prefs[3,2]))

####
##### Wrangling the GFF3 file
# ### function to pars atributes from package "davidTiling" (not the data.table way.... == slow)
# getAttributeField <- function (x, field, attrsep = ";") {
#   s = strsplit(x, split = attrsep, fixed = TRUE)
#   sapply(s, function(atts) {
#     a = strsplit(atts, split = "=", fixed = TRUE)
#     m = match(field, sapply(a, "[", 1))
#     if (!is.na(m)) {
#       rv = a[[m]][2]
#     }
#     else {
#       rv = as.character(NA)
#     }
#     return(rv)
#   })
# }
# 
# 
# gff3 <-read.gff(file="Software/Ha412HOv2.0-20181130.gff3")
# 
# gff3<-data.table(gff3)
# 
# #### now to strip the atributes
# gff3$Name <- getAttributeField(gff3$attributes, "Name")
# gff3$ID <- getAttributeField(gff3$attributes, "ID")
# gff3$Parent <- getAttributeField(gff3$attributes, "Parent")
# gff3$locus_tag <- getAttributeField(gff3$attributes, "locus_tag")
# gff3$product <- getAttributeField(gff3$attributes, "product")
# gff3$Ontology_term <- getAttributeField(gff3$attributes, "Ontology_term")
# gff3$est_cons <- getAttributeField(gff3$attributes, "est_cons")
# gff3$est_incons <- getAttributeField(gff3$attributes, "est_incons")
# gff3$ec_number <- getAttributeField(gff3$attributes, "ec_number")
# gff3$Dbxref <- getAttributeField(gff3$attributes, "Dbxref")
# 
# write.table(gff3, "Software/gff3_412HO_formated.txt",sep="\t", row.names=F, col.names=T)


gff3<-fread("Software/gff3_412HO_formated.txt")
mrna<-gff3[gff3$type=="mRNA",]
mrna$product<-url_decode(mrna$product)
mrna$attributes<-NULL



# sig.blocks<-read.table("Tables/Blocks/traits_to_genomeblocks_signif.txt", header=T)
# sug.blocks<-read.table("Tables/Blocks/traits_to_genomeblocks_sugest.txt", header=T)
### name haplotype blocks and list snps per haplotype block
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
    
    big.list<-rbind(big.list,missing.snps)

#######
sig.list<-read.table("Tables/Blocks/sigsnips_to_genomeblocks.txt",header=T)
genemap<-read.table("Tables/Blocks/condensed_genome_blocks.txt",header=T)
colocate<-read.table("Tables/Blocks/colocate_table.txt")
genemap$colocate.block<-genemap$colocate.region #rename
sig.snips<-read.table("Tables/Blocks/signif_snps_alltraits.txt")


### add colocate block name to block key
# genemap$colocate.block<-colocate$region[match(genemap$sig.hap,colocate$sighap)]


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

#### find genes in the genome haplotype blocks

genemap$nr.genes<-NA

gene.list<-NULL

for (i in 1:length(genemap$genome.hap)) {
  print(paste("Chromosome:",genemap$chr[i],"Block ID:",genemap$genome.hap[i])) 
  chrom.mrna<-mrna[mrna$seqid==genemap$chr[i]]
  
  if(sum(chrom.mrna$start>genemap$start[i]&chrom.mrna$end<genemap$stop[i])>0) {
  
      genemap$nr.genes[i]<-sum(chrom.mrna$start>genemap$start[i]&chrom.mrna$end<genemap$stop[i])
      block.mrna<-chrom.mrna[which((chrom.mrna$start>genemap$start[i]&chrom.mrna$end<genemap$stop[i])), ]
  }
  if(sum(chrom.mrna$start>genemap$start[i]&chrom.mrna$end<genemap$stop[i])==0) {
      single.snp.genes<-unique(c(last(which(chrom.mrna$start<genemap$start[i])),
                                  first(which(chrom.mrna$end>genemap$stop[i]))))
      genemap$nr.genes[i]<-length(single.snp.genes)
      block.mrna<-chrom.mrna[single.snp.genes, ]
  }
  block.mrna<-block.mrna[, c(4,5,12:18)] #remove info
  block.mrna$genome.hap<-genemap$genome.hap[i] # add genome.hap id
  block.mrna$chr<-genemap$chr[i]
  
  gene.list<-rbind(gene.list,block.mrna)  # merge together
  
}

#### wrangle with the significant snps haplotype blocks

gene.list$sig.hap<-genemap$sig.hap[match(gene.list$genome.hap,genemap$genome.hap)]
gene.list$colocate.block<-genemap$colocate.block[match(gene.list$genome.hap,genemap$genome.hap)]

gene.list<-gene.list %>% group_by(colocate.block) %>% group_by (locus_tag) %>% slice(1) ## remove duplicate genes from singif snps block

gene.list<-gene.list[,c(13,1:12)] #shuffle columns for saving


##### add traits for which the block is significant



write.csv(gene.list,"Tables/Genes/genelist.csv")

### plot region gene sizes

gene.count<-gene.list %>% group_by(colocate.block,chr) %>% count(colocate.block)

colours<-as.character(read.csv(file="Software/20colors.csv")[,2])[1:20]

gene.count$CHR<-gsub("Ha412HOChr","",gene.count$chr)

gene.count$colocate.block<-fct_rev(fct_reorder(gene.count$colocate.block, gene.count$n))

plotbase<-ggplot(gene.count,aes(x=colocate.block ,y=n, fill=CHR))
genes.plot<-plotbase+geom_point(shape=21,col="gray",size=2)+
  scale_y_continuous(trans='log10',breaks=c(1,2,3,4,5,10,50,100,500,1000))+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5,hjust = 1))+
  theme(legend.position="right")+guides(fill = guide_legend(ncol=1))+
  scale_fill_manual(values=colours,name="Chromosome")+
  xlab("significant region")+ylab("Number of genes (log10)")

genes.plot

ggsave("Plots/Colocalization/nr_genes.pdf",genes.plot, width=18, height=10)
write.csv(gene.count,"Tables/Genes/genecount.csv")
