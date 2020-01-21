##### Sunflower SAM population GWAS using GEMMA -
#####   a reimagined pipeline borowing some ideas from the version in Masalia et al 2018 (PLOS).
##### Genome: HA412HO
##### SNPs: XRQv1 SNP set built by Greg Owens (UBC) but reordered to HA412HO
##### Algorithm: GEMMA
library(data.table)

#### read in preferences
prefs<-read.table("Scripts/### Preferences ###",header=F,sep="=",skip=1)
  SNPset<-as.character(prefs[2,2])
  pheno.name<-as.character(prefs[1,2])
  multcomp<-as.numeric(as.character(prefs[3,2]))


## Read in traits and environments to run
envs<-as.character(read.table("environments_to_run.txt")[,1])
traits<-as.character(read.table("traits_to_run.txt")[,1])

pheno.data<-fread(paste("Phenotype data/",pheno.name,sep=""))


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

fam.file<-fread(paste(SNPset,".fam",sep=""))
fam.file$V6<-NULL
fam.file <- merge(fam.file,trait.data,by.x="V1",by.y="SAM",all.x=T)
write.table(file=paste(SNPset,".fam",sep=""),fam.file,col.names=F, row.names=F, quote =F)


system(paste("./gemma -bfile ",SNPset," -k ",SNPset,".cXX.txt -c ",SNPset,".PCA_EV -lmm 1 -outdir ../Tables/Assoc_files/ -o " ,paste(trait,env,sep="_"),sep=""))

  }
}

setwd("..")
