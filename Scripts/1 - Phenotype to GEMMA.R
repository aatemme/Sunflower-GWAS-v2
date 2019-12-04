##### Sunflower SAM population GWAS using GEMMA -
#####   a reimagined pipeline borowing some ideas from the version in Masalia et al 2018 (PLOS).
##### Genome: HA412HO
##### SNPs: XRQv1 SNP set built by Greg Owens (UBC) but reordered to HA412HO
##### Algorithm: GEMMA
library(data.table)

# ##### first run, convert the .tped and .tfam files into .bed, .bim, and .fam
#   wd<-getwd()
#   setwd("Software") #poor practices but typing out the path for the "system" command is a bother
#   system("./plink --tped XRQv1_412_239_filtered.tped --tfam XRQv1_412_239_filtered.tfam --recode --make-bed --out XRQv1_412_239_filtered --allow-extra-chr")
#   ### Then create the relatedness matrix
#   fam<-fread("XRQv1_412_239_filtered.fam")
#   fam$V6<-c(rep(1,dim(fam)[1]))   #### this is crucial as GEMMA will not compute relatedness if the phenotype is missing (-9)
#   write.table(file="XRQv1_412_239_filtered.fam",fam,col.names=F, row.names=F, quote =F)
#   system("./gemma -bfile XRQv1_412_239_filtered -gk 1 -o XRQv1_412_239_filtered")
#   structure<-fread("XRQv1_412_239_filtered.PCA_EV")
#   write.table(file="XRQv1_412_filtered.PCA_EV",structure ,row.names=F,col.names=F)
#   setwd(wd)
# #########


## Read in traits and environments to run
envs<-as.character(read.table("environments_to_run.txt")[,1])
traits<-as.character(read.table("traits_to_run.txt")[,1])

pheno.data<-fread("Phenotype data/GEMMAsalt.csv")
fam.file<-fread("Software/XRQv1_412_239_filtered.fam")
fam.file$V6<-NULL

setwd("Software")
for (i in 1:length(envs)){
  
  env<-envs[i]
  
  for (q in 1:length(traits)) {
    
    trait<-traits[q]
    print(paste(trait,env,sep="_"))

select_cols = c("SAM", paste(trait,env,sep="_"))

if (!select_cols[2]%in%names(pheno.data)) { 
  print("phenotype missing")
  next } #if the phenotype does not exist go to the next one

trait.data<-pheno.data[, ..select_cols]

fam.file<-fread("XRQv1_412_239_filtered.fam")
fam.file$V6<-NULL
fam.file <- merge(fam.file,trait.data,by.x="V1",by.y="SAM")
write.table(file="XRQv1_412_239_filtered.fam",fam.file,col.names=F, row.names=F, quote =F)


system(paste("./gemma -bfile XRQv1_412_239_filtered -k XRQv1_412_239_filtered.cXX.txt -c XRQv1_412_239_filtered.PCA_EV -lmm 1 -outdir ../Tables/Assoc_files/ -o " ,paste(trait,env,sep="_"),sep=""))

  }
}

setwd("..")
