#### dendrogram of traits in GWAS
library(corrr)
library(Hmisc)
library(ggdendro)
library(stringi)

##### data in the colocate plot

traits<-unique(plot.data$trait)

correlate.data<-read.table(paste("Phenotype data/",traits[1],"_",unique(plot.data$env),".txt",sep=""),header=T)
names(correlate.data)[2]<-traits[1]

for(q in 2:length(traits)){
  merge.data<-read.table(paste("Phenotype data/",traits[q],"_",unique(plot.data$env),".txt",sep=""),header=T)
  names(merge.data)[2]<-traits[q]
  
  correlate.data<-merge(correlate.data,merge.data,by="Line")
  rm(merge.data)
  
}

correlate.data$Line<-NULL

### make environment panel correlations
Env.corr<-rcorr(as.matrix(correlate.data),type="pearson")

#### clustering dendrogram to use with collocalization plot
Env.dist<-as.dist(1-Env.corr$r)    
  Env.clust <- hclust(Env.dist, method = "complete", members=NULL)
    Env.dendro<-ggdendrogram(Env.clust,rotate=T,labels=F)+scale_y_reverse()+theme(axis.text.y=element_blank())+scale_x_continuous(limits=c(0.5,length(traits)+0.5),expand=c(0,0))
      Env.label.order<-Env.clust$labels[Env.clust$order]

      
  