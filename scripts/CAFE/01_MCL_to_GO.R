setwd("~/bigdata/Suillus_comp_genomics/All_GOterms")
options(stringsAsFactors = FALSE)

#laod libraries
#library(dplyr)

#read in files
contracted_genes_in_ea_fam<- read.table("contracted_genes_in_ea_fam.txt", header = FALSE, fill = TRUE, sep = "\t")
expanded_genes_in_ea_fam<- read.table("expanded_genes_in_ea_fam.txt", header = FALSE, fill = TRUE, sep = "\t")

#read in all GO files 
#read in files
filesnames<- list.files(path ="~/bigdata/Suillus_comp_genomics/All_GOterms/", pattern = "*GO.tab", full.names=T)

#GO_files<-lapply(filesnames, fread, sep = "\t", header = TRUE, colClasses = 'character', stringsAsFactors=FALSE)

GO_files <- rbindlist(sapply(filesnames, fread, simplify = FALSE),
                use.names = TRUE, idcol = "FileName", fill = TRUE)


###format GO files
#clean up sp.names 
GO_files$species<- gsub(".*//", "", GO_files$FileName)
GO_files$species2<- gsub("_GeneCatalog.*", "", GO_files$species)

#connect sp.names to protein names
GO_files$ID <- paste(GO_files$species2, GO_files$`#proteinId`, sep="_")

#subset to only molecular function 
GO_files<- GO_files[GO_files$gotermType == "molecular_function",]

#subset to only the relevent cols (just because this is really large and it's taking a long time to process)
GO_files<- GO_files[ , c("gotermId","goName", "goAcc", "ID")]


###format expanded/contracted files
expanded_genes_in_ea_fam_formatted <- data.frame(lapply(expanded_genes_in_ea_fam, function(x) {
                 sub(".*\\|", "", x)
             }))
contracted_genes_in_ea_fam_formatted <- data.frame(lapply(contracted_genes_in_ea_fam, function(x) {
  sub(".*\\|", "", x)
}))


###get GO terms
#for each row in expanded_gnes_in_ea_fam_formatted
#if the col matches the GO_files$ID, get the GO_files$goName

#GO_df<- as.list()
#for (i in 1:nrow(expanded_genes_in_ea_fam_formatted)){
#  for (j in 1:nrow(GO_files))
#     if(expanded_genes_in_ea_fam_formatted[i] == GO_files$ID[j])
#     GO_df[i]<- as.list(GO_files$goName[j])}

#test<- setDT(expanded_genes_in_ea_fam_formatted)[V1 %chin% GO_files$ID]



###print and finish in bash 
#first replace spaces for printing in the GO file's go terms 

GO_files$goName<- GO_files[, gsub(" ", "_", GO_files$goName)]


write.table(expanded_genes_in_ea_fam_formatted, "genes_in_expanded_fams.tab", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(contracted_genes_in_ea_fam_formatted, "genes_in_contracted_fams.tab", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(GO_files, "GO_files.tab", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
