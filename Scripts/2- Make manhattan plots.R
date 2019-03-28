##### load libraries for drawing pretty manhattan plots
library(qqman)
library(tidyverse)
library(data.table)
# library(wesanderson)


###setup the data

envs<-as.character(read.table("Environments.list")[,1])
traits<-as.character(read.table("Traits.list")[,1])

suggthresh<-0.001 ## draw line at "suggestive" SNPs (threshold fraction of snips are above the blue line)

### easy map to allign the snips to
snp.map = read.table(file="Software/XRQ_264.map",header= T) ###for easier plotting of snip locations
colnames(snp.map)<-c("SNP_ID","Chr","Pos")

#### loop over every trait and then over every environment to plot the manhattans
#### loops within loops (insert Dune reference here)

i<-1
q<-1

for (i in 1:length(traits)){
  
  pdf(paste("Plots/Manhattans/",traits[i],"_ManhattanPlot.pdf",sep=""),height=8.5,width=11)
  
  if(length(envs)<3){par(mfcol=c(ceiling(length(envs)/2),1),oma = c(0, 0, 0, 0))}
              else  {par(mfcol=c(ceiling(length(envs)/2),2),oma = c(0, 0, 0, 0))} ## plot in two columns if more than 2 environments
  
  for (q in 1:length(envs)){
    
      
      snips<-fread(paste("Tables/PSfiles/",paste(traits[i],envs[q],sep="_"),".ps",sep=""),header=F)
      
      ##### bit that makes the plot
      pvalue <- snips$V4
      ytop = -log10(min(pvalue))
      if (ytop > 10){ ytop = ytop+2} else { ytop = 10}
      
      ps1<- merge(snp.map,snips[,c(1,4)],by.x=1,by.y=1)
      
      tmpcutoff <- as.data.frame(quantile(ps1$V4,as.numeric(as.character(suggthresh))))[1,1]
      
      label<-paste(traits[i],envs[q])
      print(label)

      if (length(ps1[ps1$V4<0.01,]$V4)<1) {plot(10:1)  ## draw empty plot if no snips are above minimum plotting threshold, plot fails otherwise
      }else { 
      manhattan(ps1[ps1$V4<0.01,],chr = "Chr", bp = "Pos", p = "V4", snp = "SNP_ID",col=rev(c("#3B9AB2","#E1AF00")),
                suggestiveline = -log10(tmpcutoff),genomewideline=-log10(0.05/92747), cex=0.5, cex.axis=0.8,ylim=c(2,ytop),
                mtext(label, side = 3, line = 0))
            }
      }
    mtext(label, outer = TRUE, cex = 1)
    dev.off()
      # }

}

