### single trait manhattans
##### load libraries for drawing pretty manhattan plots
library(qqman)
library(tidyverse)
library(data.table)
library(RColorBrewer)

#### read in preferences
prefs<-read.table("Scripts/### Preferences ###",header=F,sep="=",skip=1)
  SNPset<-as.character(prefs[2,2])
  pheno.name<-as.character(prefs[1,2])
  multcomp<-as.numeric(as.character(prefs[3,2]))

###setup the data

envs<-as.character(read.table("environments_to_run.txt")[,1])
traits<-as.character(read.table("traits_to_run.txt")[,1])


suggthresh<-0.001 ## draw line at "suggestive" SNPs (threshold fraction of snips are above the blue line)

### easy map to allign the snips to
snp.map = fread(file=paste("Software/",SNPset,".map",sep=""),header= F) ###for easier plotting of snip locations
colnames(snp.map)<-c("chr","rs","Chr_num","ps")
snp.map$Chr_num <- as.integer(gsub("Ha412HOChr","",snp.map$chr))

#### loop over every trait and then over every environment to plot the manhattans
#### loops within loops (insert Dune reference here)

i<-1
q<-1

for (i in 1:length(traits)){
  for (q in 1:length(envs)){
    
  pdf(paste("Plots/Manhattans/single_env/",traits[i],"-",envs[q],"_ManhattanPlot.pdf",sep=""),height=5.5,width=7.5)
  
    snips<-fread(paste("Tables/Assoc_files/",paste(traits[i],envs[q],sep="_"),".assoc.txt",sep=""),header=T)
    
    ##### bit that makes the plot
    pvalue <- snips$p_wald
    ytop <- -log10(min(pvalue))
    if (ytop > 10){ ytop = ytop+2} else { ytop = 10}
    
    ps1<- merge(snp.map,snips,by=c("chr","rs","ps"))
    
    tmpcutoff <- as.data.frame(quantile(ps1$p_wald,as.numeric(as.character(suggthresh))))[1,1]
    
    label<-paste(traits[i],envs[q])
    print(label)
    
    if (length(ps1[ps1$p_wald<0.1,]$p_wald)<1) {plot(10:1)  ## draw empty plot if no snips are above minimum plotting threshold, plot fails otherwise
    }else { 
      manhattan(ps1[ps1$p_wald<0.1,],chr = "Chr_num", bp = "ps", p = "p_wald", snp = "rs",col=brewer.pal(8,"Dark2")[c(1,8)],
                suggestiveline = -log10(tmpcutoff),genomewideline=-log10(0.05/(multcomp)), cex=0.5, cex.axis=0.8,ylim=c(1,ytop),
                mtext(label, side = 3, line = 0))
    }
  
  mtext(label, outer = TRUE, cex = 1)
  dev.off()
   }
  
}

