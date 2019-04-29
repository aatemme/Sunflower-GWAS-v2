#### dendrogram of traits in GWAS
library(corrr)
library(Hmisc)
library(ggdendro)
# library(stringi)

pheno.data<-read.csv("Phenotype data/GEMMAsalt.csv")

##### data in the colocate plot

traits<-unique(plot.data$trait)
traits.env<-paste(traits,"_",envs[i],sep="")


correlate.data<-pheno.data[,match(traits.env,names(pheno.data))]

names(correlate.data)<-traits




### make environment panel correlations
Env.corr<-rcorr(as.matrix(correlate.data),type="pearson")

#### clustering dendrogram to use with collocalization plot
Env.dist<-as.dist(1-Env.corr$r)    
  Env.clust <- hclust(Env.dist, method = "complete", members=NULL)
    Env.dendro<-ggdendrogram(Env.clust,rotate=T,labels=F)+scale_y_reverse()+theme(axis.text.y=element_blank())+scale_x_continuous(limits=c(0.5,length(traits)+0.5),expand=c(0,0))
      Env.label.order<-Env.clust$labels[Env.clust$order]

      
  