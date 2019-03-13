### This script runs the list of snips that are significant or suggestive at least once across all trait+env combinations through
### the programme LDselect per chromosome. This creates a bin of SNPs that are in LD with each other. Only bins that have a SNP that
### has been significant at least once across all trait+env combinations are retained. The output is a list of bins with their start
### and stop position and a list of SNP's that are contained in each bin. Code based on GWAS pipeline from Masalia et al (2018 PLOS)
rm(list=ls())
library(data.table)

######## LD threshold used in defining bins
ldr2 <- 0.8
########



source("Scripts/3b- Combine snips over traits and envs.R")  ## this script takes all .ps files from the EMMAX output and selects the SNPs


bin.snips<-all.snips ### rename

bin.snips$trait<-NULL
bin.snips$env<-NULL


sig.snips$Chr <- gsub("HanXRQChr","",sig.snips$V1)
sig.snips$Chr <- as.integer(gsub(":.*","",sig.snips$Chr))
sig.snips$Pos <- gsub(".*:","",sig.snips$V1)

bin.snips$Chr <- gsub("HanXRQChr","",bin.snips$V1)
bin.snips$Chr <- as.integer(gsub(":.*","",bin.snips$Chr))
bin.snips$Pos <- gsub(".*:","",bin.snips$V1)


# setwd("/home/atemme/Desktop/DREC17/EnchiladaSuite/SALSA/") ##### Change PATH



############### Loop through chromosomes present
  picmat <- NULL #For the PIC matrix at the end ## delete?
  tped <- fread(file="Software/XRQ_264.tped", header = F)
  colnames(tped)[2] = "markers"
  BIG.bins<-NULL
  
for (i in 1:length(unique(bin.snips$Chr))) {
  
  
  
  c<-sort(unique(bin.snips$Chr))[i] #select chromosome to run
  
  print(paste("Starting chromsome-",c," ",Sys.time(),sep=""))
  
  MarkMat <- data.frame()  ## delete?
  chr.snips <- bin.snips[which(bin.snips$Chr==c),]
  
  
  
  sig.chr.snips<-sig.snips[which(sig.snips$Chr==c),]
  
  print(paste(length(chr.snips$V1),"markers/", length(sig.chr.snips$V1)," significant",Sys.time()))
  
  ############################### move to next loop iteration if there are no significant snips on that chromosome
  if(length(sig.chr.snips$V1)==0) { print ("no signsnips to select" ) 
                                    next }
  ###############################
  
  
  chr.snips<-as.data.frame(chr.snips)
  sig.chr.snips<-as.data.frame(sig.chr.snips)
  
  chr.snips <- chr.snips[order(chr.snips$Pos),] ##interleaf Sig and Sug by Pos
  markers <- data.frame("markers"=chr.snips$V1)
  
#### bit that does ldselect
  binmat<-NULL


  tmp.tped <- merge(markers,tped,by = "markers")
  write.table(tmp.tped,file="Software/tmp_ldselect.tped",row.names = F,quote = F, col.names = F, sep="\t")
  
  system ("cd Software && perl ConvertToPB.pl")
  system (paste("cd Software && perl ldSelect.pl -pb Ready_tmp_ldselect.tped -freq 0.05 -r2 ",ldr2," >tmp.bins"))
  
  ldselectthing <- readLines("Software/tmp.bins") #read back in LD select results

#### bit that parses ldselects output
  #### Number of Bins including suggestive peaks #####
  NoBins = gsub("\tother_snps:.*","",ldselectthing[length(ldselectthing)-1])
  NoBins = gsub("Bin ","",NoBins)
  NoBins = as.numeric(as.character(gsub(" ","",NoBins)))
  
  marker_bin_mat = NULL
  tmp_bin_data = NULL
  SSBinTable = NULL
  ### create table of bins and markers ####
  sites = grep("total_sites:",readLines("Software/tmp.bins"))
  
  
  for (i in 1:length(sites)){
    tmp_site = sites[i]
    tmpbin = gsub("\t.*","",ldselectthing[tmp_site])
    tmpbin = as.numeric(as.character(gsub("Bin ","",tmpbin)))
    tmp_nosites = gsub(".*\ttotal_sites: ","",ldselectthing[tmp_site])
    tmp_nosites = as.numeric(as.character(gsub("\t.*","",tmp_nosites)))
    tmp_row = cbind(tmpbin,tmp_nosites)
    binmat = rbind(binmat,tmp_row)
    
    tagline = tmp_site+1
    tmptag = gsub(".*: ","",ldselectthing[tagline])
    tmptag = gsub("\\s+$","",tmptag,perl=T)
    tmptag = strsplit(tmptag," ",perl=T) 
    tmptag = unlist(tmptag)
    if(length(tmptag)==0){next}
    else {
      for(b in 1:length(tmptag)){
        blah = cbind(as.numeric(as.character(tmptag[b])),tmpbin)
        marker_bin_mat = rbind(marker_bin_mat,blah)
      }
    }
    osline = tmp_site+2
    tmpos = gsub(".*: ","",ldselectthing[osline])
    tmpos = gsub("\\s+$","",tmpos,perl=T)
    tmpos = strsplit(tmpos," ",perl=T) 
    tmpos = unlist(tmpos)
    if(length(tmpos)==0){next}
    else {
      for(s in 1:length(tmpos)){
        blahos = cbind(as.numeric(as.character(tmpos[s])),tmpbin)
        marker_bin_mat = rbind(marker_bin_mat,blahos)
      }
      binblah = cbind(tmpbin,as.numeric(as.character(tmpos[1])),as.numeric(as.character(tmpos[length(tmpos)])))
      SSBinTable = rbind(SSBinTable,binblah)
    }
  } #End length sites
  marker_bin_mat<-as.data.frame(marker_bin_mat)
  colnames(marker_bin_mat) = c("Pos","Bin")

####  
  
  ### write some bin outputs ###
  colnames(binmat) = c("Bin","No SNPs Per Bin")
  #write.table(binmat,file="NoBins_SNPsPerBin.list",row.names = F,quote = F, col.names = T, sep="\t")
  
  ### Rematch with Sig markers:
  tmp_bin_data <- merge(sig.chr.snips, marker_bin_mat, by="Pos")
  tmp_bin_data = tmp_bin_data[order(tmp_bin_data$Pos),]
  BinsofInterest = levels(as.factor(tmp_bin_data$Bin))
  
  tmp_chrom_data = merge(chr.snips,marker_bin_mat,by="Pos")
  tmp_chrom_data = tmp_chrom_data[order(tmp_chrom_data$Pos),]
  
  RegionsofInterest = tmp_chrom_data[tmp_chrom_data$Bin%in%BinsofInterest, ]
  

  RegionsofInterest$SNP_ID<-RegionsofInterest$V1
  
  #backup = RegionsofInterest
  print(paste("ldselect done ",length(BinsofInterest)," regions",Sys.time()))
###### see if bin is inbetween another bin #####
  test = RegionsofInterest[,c("SNP_ID","Bin")]
  
  test$Bin  = as.factor(test$Bin)
  test$pos = as.numeric(as.character((gsub(".*:","",test$SNP_ID))))
  test<-test[order(test$pos),]
  
  binss = NULL
  for(o in 1:length(levels(test$Bin))){
    jose = levels(test$Bin)[o]
    subbin = test[which(test$Bin==jose),]
    binstart = subbin$pos[1]
    binstop = subbin$pos[length(subbin$pos)]
    bdiff = binstop - binstart
    btmp = cbind(o,binstart,binstop,bdiff)
    binss = rbind(binss,btmp)
  }
  binss = as.data.frame(binss)
  binss = binss[order(rev(binss$bdiff)),]
  sharpie=NULL
  for(y in 1:nrow(test)){
    checkmark = test$pos[y]
    pen = NULL
    for (z in 1:nrow(binss)){
      bcstart = binss[z,2]
      bcstop =  binss[z,3]
      if (checkmark >= bcstart && checkmark <= bcstop){checkbin = binss[z,1]
      } else {checkbin = test[y,3]}
      checktmp = cbind(checkmark,checkbin)
      pen = rbind(pen,checktmp)
      pen = as.data.frame(pen)
      trubin = min(pen$checkbin)
      pentmp = cbind(checkmark,trubin)
    }
    sharpie = rbind(sharpie,pentmp)
  }
  
  test = as.data.frame(sharpie)
  colnames(test)[1] = "Pos"
  test = merge(test,RegionsofInterest, by ="Pos")  
  
  test$pretty.bin<-as.numeric(as.factor(test$trubin))
  test<-test[order(test$Pos), ]
  
  print(paste("Regions colapsed to ",max(test$pretty.bin),Sys.time()))


  BIG.bins<-rbind(BIG.bins,test) ### grow list of snips in bins with each iteration of the chromosome loop
  
  write.csv(x=test,paste("Tables/Bins/chromosome-",c,".csv",sep="")) ### store bins per chromosome
 
print(paste("chromsome-",c," done ",Sys.time(),sep=""))  

### cleanup
system("rm Software/tmp.bins")
system("rm Software/Ready_tmp_ldselect.tped")
system("rm Software/tmp_ldselect.tped")
   
}  
  
save<-BIG.bins  
  
BIG.bins$region<-as.factor(paste(BIG.bins$Chr,BIG.bins$pretty.bin,sep="_")) 

write.csv(x=BIG.bins,paste("Tables/Bins/All-bins.csv",sep=""))


picmat<- NULL  
  
  for(p in 1:length(levels(BIG.bins$region))){
    mouse = as.character(levels(BIG.bins$region)[p])
    subset = BIG.bins[which(BIG.bins$region==mouse),]
    start = as.character(subset[1,1])
    stop = as.character(subset[nrow(subset),1])
    ds = "across all"
    chromosome = subset$Chr[1]
    binnybinbin <- mouse
    picline = cbind(ds,chromosome,start,stop, binnybinbin)
    
    picmat = rbind(picmat,picline)
  }

write.csv(x=picmat,paste("Tables/Bins/All-regions.csv",sep=""))

picmat<-as.data.frame(picmat)
picmat$chromosome<-as.character(picmat$chromosome)

picmat$chrom<-formatC(as.numeric(picmat$chromosome),width=2, flag="0")
picmat$start<-paste("HanXRQChr",picmat$chrom,":",picmat$start,sep="")
picmat$stop<-paste("HanXRQChr",picmat$chrom,":",picmat$stop,sep="")
picmat$chrom<-NULL


### file for creating the gene list
write.table(picmat,file="Tables/Bins/Pre_Identified_Clusters.csv",row.names = F,quote = F, col.names = T, sep="\t")



    