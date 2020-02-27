setwd("~/bigdata/Suillus_comp_genomics/All_GOterms/fams_of_interest")
options(stringsAsFactors = FALSE)

#laod libraries
library(dplyr)
library(data.table)

#read in files
contracted_genes_in_ea_fam<- read.table("contracted_genes_in_ea_fam.txtt", header = FALSE, fill = TRUE, sep = "\t")
expanded_genes_in_ea_fam<- read.table("expanded_genes_in_ea_fam.txtt", header = FALSE, fill = TRUE, sep = "\t")

#format expanded/contracted files
expanded_genes_in_ea_fam_formatted <- data.frame(lapply(expanded_genes_in_ea_fam, function(x) {
  sub(".*\\|", "", x)
}))
contracted_genes_in_ea_fam_formatted <- data.frame(lapply(contracted_genes_in_ea_fam, function(x) {
  sub(".*\\|", "", x)
}))

#read in the Interpro terms files 
#read in files
filesnames<- list.files(path ="~/bigdata/Suillus_comp_genomics/All_InterPro/", pattern = "*IPR.tab", full.names=T)

#rename the input speices path to just species
Interpro_files <- rbindlist(sapply(filesnames, fread, simplify = FALSE),
                            use.names = TRUE, idcol = "FileName", fill = TRUE)
Interpro_files2 <- data.frame(lapply(Interpro_files, function(x) {
  sub(".*//|", "", x)
}))
Interpro_files3 <- data.frame(lapply(Interpro_files2, function(x) {
  sub("_Gene.*", "", x)
}))

#remove lines with no annotations
Interpro_files3<- Interpro_files3[! Interpro_files3$iprId == "\\N",]

#combine Sp name and prot ID cols into a new col. 
Interpro_files3$ID<- paste(Interpro_files3$FileName, Interpro_files3$X.proteinId, sep="_")

#shrink to just the cols of interest
IPR_of_interest<- Interpro_files3[ , c("iprId","iprDesc", "domainDb", "ID")]

#reformat spaces for printing 
IPR_of_interest$iprDesc<- gsub(",", "|", IPR_of_interest$iprDesc)
IPR_of_interest$iprDesc<- gsub("\\| ", "|", IPR_of_interest$iprDesc)
IPR_of_interest$iprDesc<- gsub(" ", "_", IPR_of_interest$iprDesc)


write.table(IPR_of_interest, "IPR_files.tab", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)



