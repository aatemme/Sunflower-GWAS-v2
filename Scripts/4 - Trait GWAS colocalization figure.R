library(gridExtra)
library(ggpubr)
library(cowplot)
library(wesanderson)

colors<-c(wes_palette("Darjeeling1")[1],wes_palette("Darjeeling1")[5])  

envs<-as.character(read.table("environments_to_run.txt")[,1])

sig.blocks<-read.table("Tables/Blocks/traits_to_genomeblocks_signif.txt", header=T)
sug.blocks<-read.table("Tables/Blocks/traits_to_genomeblocks_sugest.txt", header=T)
sig.list<-read.table("Tables/Blocks/sigsnips_to_genomeblocks.txt",header=T)
sighap_to_genomehap<-read.table("Tables/Blocks/condensed_genome_blocks.txt",header=T)

colocate<-rbind(sig.blocks,sug.blocks[sug.blocks$hapID%in%sig.blocks$hapID,])

colocate$sighap<-sighap_to_genomehap$sig.hap[match(colocate$hapID,sighap_to_genomehap$genome.hap)]

colocate$region<-sighap_to_genomehap$colocate.region[match(colocate$hapID,sighap_to_genomehap$genome.hap)]

colocate$trait_env<-paste(colocate$trait,colocate$env,sep="_")

traits.per.block<-colocate %>% group_by(sighap) %>% summarise(trait_num=length(trait_env))

# single.trait.blocks<-colocate[colocate$hapID%in%traits.per.block$hapID[traits.per.block$trait_num==1],]
# 
# colocate<-colocate[!colocate$hapID%in%traits.per.block$hapID[traits.per.block$trait_num==1],]

colocate<-colocate %>% separate(sighap, sep= "_", c("chromosome","blocknum"),remove=F) %>%
                        arrange(chromosome, blocknum)

colocate<- colocate %>% mutate(beta.sign=sign(beta))

colocate$region<-factor(colocate$region)

write.table(colocate,"Tables/Blocks/colocate_table.txt")

chrom.borders<-colocate %>% group_by(chromosome)%>% summarise(bin=length(unique(region))) %>% arrange(as.integer(chromosome))
chrom.borders<-cumsum(chrom.borders$bin)
chrom.borders<-chrom.borders+0.5
chrom.borders<-chrom.borders[1:length(chrom.borders)-1]

# baseplot<-ggplot(colocate,aes(x=region,y=trait,fill=as.factor(beta.sign)))
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
#   ggtitle("control")+theme(legend.position = "none")+theme(axis.title.x=element_blank())+
#   facet_wrap(~env,nrow=3)

#### draw the tree environment colocate plots separately

for (i in 1: length(envs)) {

plot.data<-colocate[colocate$env==envs[i],]

source("Scripts/4b - correlation dendrogram.R") ### update script referal after changing names

plot.data$trait.factor<-factor(plot.data$trait,levels=Env.label.order)

baseplot<-ggplot(plot.data,aes(x=region,y=trait.factor,fill=as.factor(beta.sign)))

plot.colocate<- baseplot+geom_vline(xintercept=c(1:length(plot.data$region)),colour="darkgrey",linetype=3)+
  geom_vline(xintercept=chrom.borders,colour="black")+
  geom_tile(fill="white")+
  geom_tile(aes(alpha=pvalue),colour="black")+
  theme_minimal()+
  theme(axis.text.y = element_text(hjust = 0))+
  scale_fill_manual(values=c(colors[1],colors[2]))+
  scale_alpha_manual(values=c(1,0.1))+
  scale_x_discrete(drop=F)+
  theme_classic()+
  theme(axis.title.y=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5,hjust=1))+
  ggtitle(envs[i])+theme(legend.position = "none")+theme(axis.title.x=element_blank())

# plot.data.dendro<-plot.data.dendro+coord_flip(xlim=c(4,length(plot.data.label.order)-4))+
#   theme(axis.text.x=element_text(size=8))+theme(axis.title.x=element_blank())

comb.plot<-plot_grid(Env.dendro+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
                     plot.colocate+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
                     align="h",rel_widths=c(1,9))

trait.to.region.ratio<-length(levels(plot.data$region))/length(Env.label.order)

ggsave(paste("Plots/Colocalization/colocate-",envs[i],".pdf",sep=""),plot=comb.plot,width=12,height=floor(14/trait.to.region.ratio))

}

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



