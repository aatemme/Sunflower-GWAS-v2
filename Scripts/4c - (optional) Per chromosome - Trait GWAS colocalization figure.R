library(gridExtra)
library(ggpubr)
library(cowplot)
library(tidyverse)
library(wesanderson)

colors<-c("#1b9e77", "gray85") 

envs<-as.character(read.table("environments_to_run.txt")[,1])

sig.blocks<-read.table("Tables/Blocks/traits_to_genomeblocks_signif.txt", header=T)
sug.blocks<-read.table("Tables/Blocks/traits_to_genomeblocks_sugest.txt", header=T)
sig.list<-read.table("Tables/Blocks/sigsnips_to_genomeblocks.txt",header=T)
sighap_to_genomehap<-read.table("Tables/Blocks/condensed_genome_blocks.txt",header=T)

#### set up data to feed into plotting
colocate<-rbind(sig.blocks,sug.blocks[sug.blocks$hapID%in%sig.blocks$hapID,])
  colocate$sighap<-sighap_to_genomehap$sig.hap[match(colocate$hapID,sighap_to_genomehap$genome.hap)]
    colocate$region<-sighap_to_genomehap$colocate.region[match(colocate$hapID,sighap_to_genomehap$genome.hap)]
      colocate$trait_env<-paste(colocate$trait,colocate$env,sep="_")

  traits.per.block<-colocate %>% group_by(sighap) %>% summarise(trait_num=length(trait_env))

    colocate<-colocate %>% separate(sighap, sep= "_", c("chromosome","blocknum"),remove=F) %>%
                            arrange(chromosome, blocknum)

    colocate<- colocate %>% mutate(beta.sign=sign(beta))

    colocate$region<-factor(colocate$region)


############## function that helps with settting a region as significant if it's only significant for one of it's component genome haplotype blocks
sig.sug.fun<-function (x) {
 if (sum("significant"%in%x)>0) {y<-"signficicant"} 
  if (sum("significant"%in%x)==0) {y<-"suggestive"}
  return(y)
}
##################

##### condense to single entry per region (collapse genome blocks) 
colocate<-colocate %>% group_by(region,trait_env) %>% dplyr::summarize(trait=trait[1],  
                                                            env=env[1],
                                                            pvalue=factor(sig.sug.fun(pvalue)),
                                                            chromosome=chromosome[1],
                                                            beta.sign=factor(sign(mean(beta.sign))))

#### plot per chromosome in a loop

for (i in 1: length(unique(colocate$chromosome))) {
  
  chrom.data<-colocate[colocate$chromosome==unique(colocate$chromosome)[i],]
  chrom.data$region<-as.character(chrom.data$region)
  chrom.data$region<-factor(chrom.data$region)
  plots<-list()
  
  region.count<-length(levels(chrom.data$region))
  
  trait.count<-c(NA,NA,NA)
  
  for (q in 1:length(envs)) {
    plot.data<-chrom.data[chrom.data$env==envs[q],]
  
      source("Scripts/4b - correlation dendrogram.R") 
      
      plot.data$trait.factor<-factor(plot.data$trait,levels=Env.label.order)
      
      baseplot<-ggplot(plot.data,aes(x=region,y=trait.factor,fill=pvalue))
      
      plot.colocate<- baseplot+geom_vline(xintercept=c(1:length(plot.data$region)),colour="darkgrey",linetype=3)+
        geom_tile(fill="white")+
        geom_tile(colour="black")+
        geom_point(aes(shape=beta.sign),size=3.5)+
        scale_shape_manual(values=c("+","-","±"))+
        theme_minimal()+
        theme(axis.text.y = element_text(hjust = 0))+
        scale_fill_manual(values=c(colors[1],colors[2]))+
        scale_alpha_manual(values=c(1,0.1))+
        scale_x_discrete(drop=F)+
        theme_classic()+
        theme(axis.title.y=element_blank(),axis.text.x = element_text(angle = 90, vjust = 0.5,hjust=1))+
        ggtitle(paste(envs[q],"Chr",unique(colocate$chromosome)[i]))+theme(legend.position = "none")+theme(axis.title.x=element_blank())
      
      # plot.data.dendro<-plot.data.dendro+coord_flip(xlim=c(4,length(plot.data.label.order)-4))+
      #   theme(axis.text.x=element_text(size=8))+theme(axis.title.x=element_blank())
      
      # comb.plot<-plot_grid(Env.dendro+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
      #                      plot.colocate+theme(plot.margin = unit(c(0, 0, 0, 0), "cm")),
      #                      align="h",rel_widths=c(1,9))
      
      Env.dendro<-Env.dendro+theme(axis.text.x = element_text(angle = 90, vjust = 0.5,hjust=1,size=8))+ggtitle(" ")+theme(plot.margin = unit(c(0.2, 0, 0, 0), "cm"))
      
      assign(envs[q],plot.colocate)
      assign(paste(envs[q],"_dendro",sep=""),Env.dendro)
      trait.count[q]<-length(Env.label.order)
  }
  
  trait.count<-trait.count[c(2,1,3)]
  trait.plot<-plot_grid(water,salt,logdiff,align="v",nrow=3,rel_heights=trait.count+2.5)
  dendro.plot<-plot_grid(water_dendro,salt_dendro,logdiff_dendro,align="v",nrow=3,rel_heights=trait.count+2.5)
  
  chrom.plot<-plot_grid(dendro.plot,trait.plot,ncol=2,rel_widths = c(3,region.count),align="v")
  ggsave(paste("Plots/Colocalization/colocate-chromosome-",unique(colocate$chromosome)[i],".pdf",sep=""),plot=chrom.plot,width=6,height=9)
  
}

