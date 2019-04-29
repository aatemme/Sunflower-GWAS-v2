##### Wrangling the GFF3 file
library(ape)
library(data.table)


### function to pars atributes from package "davidTiling" (not the data.table way.... == slow)
getAttributeField <- function (x, field, attrsep = ";") {
  s = strsplit(x, split = attrsep, fixed = TRUE)
  sapply(s, function(atts) {
    a = strsplit(atts, split = "=", fixed = TRUE)
    m = match(field, sapply(a, "[", 1))
    if (!is.na(m)) {
      rv = a[[m]][2]
    }
    else {
      rv = as.character(NA)
    }
    return(rv)
  })
}


gff3 <-read.gff(file="Software/Ha412HOv2.0-20181130.gff3")

gff3<-data.table(gff3)

#### now to strip the atributes
gff3$Name <- getAttributeField(gff3$attributes, "Name")
gff3$ID <- getAttributeField(gff3$attributes, "ID")
gff3$Parent <- getAttributeField(gff3$attributes, "Parent")
gff3$locus_tag <- getAttributeField(gff3$attributes, "locus_tag")
gff3$product <- getAttributeField(gff3$attributes, "product")
gff3$Ontology_term <- getAttributeField(gff3$attributes, "Ontology_term")
gff3$est_cons <- getAttributeField(gff3$attributes, "est_cons")
gff3$est_incons <- getAttributeField(gff3$attributes, "est_incons")
gff3$ec_number <- getAttributeField(gff3$attributes, "ec_number")
gff3$Dbxref <- getAttributeField(gff3$attributes, "Dbxref")


mrna<-gff3[gff3$type=="mRNA",]


