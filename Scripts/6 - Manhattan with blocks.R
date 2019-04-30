### single trait manhattans with significant blocks overlay
##### load libraries for drawing pretty manhattan plots
# library(qqman)
library(tidyverse)
library(data.table)
library(RColorBrewer)
library(ggpubr)

colours<-rep(c(brewer.pal(8,"Dark2"))[-7],3)
# library(wesanderson)

sig.list<-read.table("Tables/Blocks/sigsnips_to_genomeblocks.txt",header=T)
genemap<-read.table("Tables/Blocks/condensed_genome_blocks.txt",header=T)
colocate<-read.table("Tables/Blocks/colocate_table.txt")

colocate<- colocate %>% group_by(chromosome) %>% mutate(region_col=colours[as.numeric(factor(rank(match(region,levels(region)))))])

### name haplotype blocks and list snps per haplotype block
blocks<-fread("Software/XRQv1_412_239_filtered.blocks.det")
blocks$Chr_num<- as.integer(gsub("Ha412HOChr","",blocks$CHR))
blocks<- blocks %>% group_by(Chr_num) %>% mutate(hapID = paste(Chr_num,c(1:length(Chr_num)),sep="_"))
snps<-strsplit(blocks$SNPS,split="|",fixed=T)
big.list<-unlist(snps)
big.list<-data.table(SNP=big.list)
big.list$hapID<-c(rep(blocks$hapID, blocks$NSNPS))
rm(snps)

### qd solution singletons
all.snps<-fread("Tables/Assoc_files/Calcium_logdiff.assoc.txt", header=T)

missing.snps<-all.snps[!all.snps$rs%in%big.list$SNP, c(1:3) ]
missing.snps$Chr_num<- as.integer(gsub("Ha412HOChr","",missing.snps$chr))
missing.snps<- missing.snps %>% group_by(Chr_num) %>% mutate(hapID=paste(Chr_num,"_single",match(rs,unique(rs)),sep=""))
missing.snps<-missing.snps[,c(2,5)]
names(missing.snps)<-c("SNP","hapID")

big.list<-rbind(big.list,missing.snps)

###setup the data

envs<-as.character(read.table("environments_to_run.txt")[,1])
traits<-as.character(read.table("traits_to_run.txt")[,1])


suggthresh<-0.001 ## draw line at "suggestive" SNPs (threshold fraction of snips are above the blue line)

### easy map to allign the snips to
# snp.map = fread(file="Software/XRQv1_412_239_filtered.map",header= F) ###for easier plotting of snip locations
# colnames(snp.map)<-c("chr","rs","Chr_num","ps")
# snp.map$Chr_num <- as.integer(gsub("Ha412HOChr","",snp.map$chr))

#### loop over every trait and then over every environment to plot the manhattans
#### loops within loops (insert Dune reference here)

i<-26
q<-1

# for (i in 1:length(traits)){
#   for (q in 1:length(envs)){
    
    # pdf(paste("Plots/Manhattans/single_env/",traits[i],"-",envs[q],"_ManhattanPlot.pdf",sep=""),height=5.5,width=7.5)
    label<-paste(traits[i],envs[q])
    print(label)
    snips<-fread(paste("Tables/Assoc_files/",paste(traits[i],envs[q],sep="_"),".assoc.txt",sep=""),header=T)
    snips$CHR<- as.integer(gsub("Ha412HOChr","",snips$chr))
    
    snips<-merge(snips,big.list,by.x="rs",by.y="SNP")
    
    trait.blocks<-colocate[colocate$trait_env==paste(traits[i],envs[q],sep="_"),]
    
    
    
    tmpcutoff <- -log10(as.data.frame(quantile(snips$p_wald,as.numeric(as.character(suggthresh))))[1,1])
    ##### bit that makes the plot
    
    spacer<-30000000
    
    chr_cumsum<- snips %>% ## code from https://www.r-graph-gallery.com/wp-content/uploads/2018/02/Manhattan_plot_in_R.html
            group_by(CHR) %>% 
            summarise(chr_len=max(ps)) %>%  # Compute chromosome size
            mutate(tot=cumsum(chr_len)-chr_len) %>% # Calculate cumulative position of each chromosome
            select(-chr_len)
    
    chr_cumsum$tot<-chr_cumsum$tot+c((spacer*chr_cumsum$CHR)-spacer)
    
   snips<- merge(snips, chr_cumsum, by="CHR")
   snips$BPcum<-snips$ps+snips$tot 
   
    snips<-data.table(snips)
    
    genome.size<-max(snips$BPcum)
    
    axisdf <- snips %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 ) ## set chromosome label position
    
    ytop <- -log10(min(snips$p_wald))
    if (ytop > 10){ ytop = ytop+2} else { ytop = 10}
    
    snips<-snips %>%filter(-log10(p_wald)>1)
    
    highlights<-snips[snips$hapID%in%trait.blocks$hapID,] # snips to highlight for their blocks
    
    highlights<-merge(highlights,trait.blocks,by="hapID",all.x=T)
    # highlights<-highlights %>% group_by(chromosome) %>% mutate(region_col=colours[as.numeric(factor(rank(match(region,levels(region)))))])
    highlights$region_col<-colocate$region_col[match(highlights$region,colocate$region)]
    
    
    plot<-ggplot(data=snips, aes(x=BPcum, y=-log10(p_wald),color=as.factor(CHR)))+
      geom_point(size=0.4)+scale_color_manual(values=rep(c("grey75",brewer.pal("Blues",n=9)[3]),17))+
      annotate("point",x=highlights$BPcum,y=-log10(highlights$p_wald),col=highlights$region_col,size=0.6)+
      scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center, expand = expand_scale(mult = c(0.02, 0.02))) +
      scale_y_continuous(expand = c(0, 0), limits=c(1,ytop), breaks=seq(from=2, to=ytop,by=2) )+
      theme_light() +
      theme( 
        legend.position="none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )+
      xlab("Chromosome")+ylab(expression("-log"[10] * "(p)"))+
      geom_hline(yintercept = -log10(0.05/(31733+15809)), col="red")+
      geom_hline(yintercept=tmpcutoff,col="blue")+
      ggtitle(label)
    
    ggsave(paste("Plots/Manhattans_regionhighlight/",traits[i],"_",envs[q],".png",sep=""),plot, height=4.5,width=7.5, units="in",dpi=300)
  
#   
# }
