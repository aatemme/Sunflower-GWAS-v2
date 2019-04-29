##### make synthetic trait that can be run through salsa to get the overall bins data.
library(data.table)


###setup the data

envs<-as.character(read.table("environments_to_run.txt")[,1])
traits<-as.character(read.table("traits_to_run.txt")[,1])


thresh<-0.05/(31733+15809)
suggthresh<-0.001

sig.snips<-NULL
sug.snips<-NULL


for (i in 1:length(traits)){
  for(q in 1:length(envs)) {

print(paste(traits[i],envs[q]))

snips.to.merge<-fread(paste("Tables/Assoc_files/",paste(traits[i],envs[q],sep="_"),".assoc.txt",sep=""),header=T)
    
if(range(snips.to.merge$p_wald)[1]<thresh) {
  print("sig snips!")
  tmpcutoff <- as.data.frame(quantile(snips.to.merge$p_wald,as.numeric(as.character(suggthresh))))[1,1]
  
  sig.snips.to.merge<-snips.to.merge[which(snips.to.merge$p_wald<thresh), ]
    sig.snips.to.merge$trait<-traits[i]
      sig.snips.to.merge$env<-envs[q]
  
  sug.snips.to.merge<-snips.to.merge[which(snips.to.merge$p_wald<tmpcutoff&snips.to.merge$p_wald>thresh), ]
    sug.snips.to.merge$trait<-traits[i]
      sug.snips.to.merge$env<-envs[q]
  
  
  sig.snips<-rbind(sig.snips,sig.snips.to.merge,make.row.names = FALSE)
  sug.snips<-rbind(sug.snips,sug.snips.to.merge,make.row.names = FALSE)

}
  }
}

sig.snips.doubles<-sig.snips
sug.snips.doubles<-sug.snips

sig.snips <- unique(sig.snips,by="rs")
sug.snips <- unique(sug.snips,by="rs")

all.snips<-unique(rbind(sig.snips,sug.snips),by="rs")
