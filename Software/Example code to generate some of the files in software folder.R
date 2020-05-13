#### script with some example code on how to generate some of the files needed in the software folder
####
#### NOT to be run in using the pipeline


### how to generate a kinship matrix using gemma
#   fam<-fread("XRQv1_412_239_filtered.fam")
#   fam$V6<-c(rep(1,dim(fam)[1]))   #### this is crucial as GEMMA will not compute relatedness if the phenotype is missing (-9)
#   write.table(file="XRQv1_412_239_filtered.fam",fam,col.names=F, row.names=F, quote =F)
#   system("./gemma -bfile XRQv1_412_239_filtered -gk 1 -o XRQv1_412_239_filtered")


### how to generate the haplotype blocks using plink
# plink command run per chromosome (simultaneously since it's single core)
# ./plink --tped XRQv1_412_239_filtered.tped --tfam XRQv1_412_239_filtered.tfam --blocks 'no-pheno-req' 'no-small-max-span' --blocks-max-kb 100000 --blocks-strong-lowci 0.7005 --out CHR1_infFRAC_9 --allow-extra-chr --chr Ha412HOChr01 --blocks-inform-frac 0.9
# Then combine the .blocks.det files per chromosome into one
