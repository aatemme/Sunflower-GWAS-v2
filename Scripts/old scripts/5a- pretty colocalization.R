### pretty plot of colocalization
rm(list=ls())
library(tidyverse)
library(data.table)
library(gridExtra)
library(ggpubr)
library(wesanderson)
library(cowplot)



BIG.bins<-read.csv("Tables/Bins/All-bins.csv")
  BIG.bins<-BIG.bins[,!names(BIG.bins)%in%c("V2","V3","V4","X","Pos","Bin")]
  
  
  
Regions<-read.csv("Tables/Bins/All-regions.csv")
  Regions$gene.region<-c(1:length(Regions[,1 ]))
  
  
envs<-as.character(read.table("Environments.list")[,1])
traits<-as.character(read.table("Traits.list")[,1])


  
  thresh<-0.05/92747
  suggthresh<-0.001
  
  colocate<-NULL #empty object to merge against
  

  for (i in 1:length(traits)){
    for(q in 1:length(envs)) {
      
      # if(i==56&q==4) {break} ## Janky example code to filter out a .ps file if its causing trouble
      
      print(paste(traits[i],envs[q]))
      sig.bins<-NULL
      sug.bins<-NULL
      
  snips<-fread(paste("Tables/PSfiles/", paste(traits[i],envs[q],sep="_"), ".ps", sep=""), header=F)
  tempcutoff <- as.data.frame(quantile(snips$V4, as.numeric(as.character(suggthresh))))[1, 1]
  sug.snips<-snips[which(snips$V4<tempcutoff&snips$V4>thresh), ]
  
  if (range(snips$V4)[1]<thresh){
  sig.snips<-snips[which(snips$V4<thresh), ]


  sig.bins<-merge(sig.snips,BIG.bins,by="V1")
  sig.bins<-as.data.frame(unique(sig.bins,by="region")) ### change to take median of beta
  sig.bins<-sig.bins[ ,which(names(sig.bins)%in%c("Chr","pretty.bin","region","V2"))]
  sig.bins$trait<-traits[i]
  sig.bins$env<-envs[q]
  sig.bins$pvalue<-"significant"
  }
  
 
  sug.bins<-merge(sug.snips,BIG.bins,by="V1")
  sug.bins<-as.data.frame(unique(sug.bins,by="region")) ### change to take median of beta
  sug.bins<-sug.bins[ ,which(names(sug.bins)%in%c("Chr","pretty.bin","region","V2"))]
 if(length(sug.bins$Chr)>0){
   sug.bins$trait<-traits[i]
  sug.bins$env<-envs[q]
  sug.bins$pvalue<-"suggestive"
 }
  if(length(sug.bins$Chr)==0){
    sug.bins$trait<-NULL
    sug.bins$env<-NULL
    sug.bins$pvalue<-NULL
  }
  
  
  sug.bins<-sug.bins[!sug.bins$region%in%sig.bins$region, ]
  
  
  trait.bins<-rbind(sig.bins,sug.bins)  
  
  colocate<-rbind(colocate,trait.bins)
    }
  }
 
  colocate$trait.env<-fct_rev(fct_inorder(as.factor(paste(colocate$trait,colocate$env)), ordered=T))
  colocate$pvalue<-as.factor(colocate$pvalue)
  
  colocate$region<-as.character(colocate$region)

  colocate<-colocate[
    with(colocate, order(colocate$Chr, colocate$pretty.bin)),
    ]
  
  colocate$region<-fct_inorder(colocate$region, ordered = NA)

  
  colocate$beta.sign<-as.factor(sign(colocate$V2))
  
##### save colocalization data
  write.csv(colocate,"Tables/Genes/colocate.csv")  
#######
  
###### plot colocalization per environment (note, bins are the same across all environments)  
   
  chrom.borders<-colocate %>% group_by(Chr)%>% summarise(bin=n_distinct(pretty.bin))
  chrom.borders<-cumsum(chrom.borders$bin)
  chrom.borders<-chrom.borders+0.5
  chrom.borders<-chrom.borders[1:length(chrom.borders)-1]


### split environments. Create data.frame with just the regions per environment
  for (i in 1:length(envs)) {
  
  filtered <- colocate %>%
    filter(env==envs[i]) 
  
  assign(envs[i],filtered) ## assign/get are unloved online...but it works, so...?
  
  rm(filtered)
  
  }
  


  colors<-c(wes_palette("Darjeeling1")[1],wes_palette("Darjeeling1")[5])  
  
### Draw plot for each environment
  
  for (i in 1:length(envs)) {
    
    plot.data<-get(envs[i]) ## assign/get are unloved online...but it works, so...?
    
    source("Scripts/5b- correlation dendrogram.R")
  
  plot.data$trait.factor<-factor(plot.data$trait,levels=Env.label.order)
  
  baseplot<-ggplot(plot.data,aes(x=region,y=trait.factor,fill=beta.sign))

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
  
  trait.to.region.ratio<-length(Regions[,1])/length(Env.label.order)
  
 ggsave(paste("Plots/Colocalization/colocate-",envs[i],".pdf",sep=""),plot=comb.plot,width=12,height=floor(12/trait.to.region.ratio))
  
  } 

  
  
  
  
###### plot number of genes per region  
genes<-read.table("Tables/Genes/NoGenesPerBin.txt",sep="\t",header=T)
  names(genes)<-c("gene.region","Chr","n")


colours<-as.character(read.csv(file="Software/20colors.csv")[,2])[1:20]

  region.names<-data.frame(Bin=Regions$gene.region,gene.region=Regions$binnybinbin)

  genes<-merge(genes,region.names, by="gene.region")
    genes$gene.region<-fct_rev(fct_reorder(genes$gene.region, genes$n))

  plotbase<-ggplot(genes,aes(x=gene.region,y=n, fill=as.factor(Chr)))
genes.plot<-plotbase+geom_bar(stat="identity")+
  scale_y_continuous(trans='log10',breaks=c(2,3,4,5,10,50,100,500,2000))+
  theme_minimal()+
  theme(axis.title.y=element_blank(),axis.text.x = element_text(angle = 90, vjust=0.5,hjust = 1))+
  theme(legend.position="bottom")+guides(fill = guide_legend(nrow=1))+
  scale_fill_manual(values=colours,name="Chromosome")+theme(axis.title.x=element_blank())

ggsave("Plots/Colocalization/genes number.pdf",plot=genes.plot,width=10,height=7) 


