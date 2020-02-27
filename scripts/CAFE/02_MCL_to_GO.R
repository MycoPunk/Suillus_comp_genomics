setwd("~/bigdata/Suillus_comp_genomics/All_GOterms/fams_of_interest")
options(stringsAsFactors = FALSE)

#laod libraries
#library(dplyr)

#read in files
contracted_genes_in_ea_fam<- read.table("contracted_genes_in_ea_fam.txt", header = FALSE, fill = TRUE, sep = "\t", row.names=1)
expanded_genes_in_ea_fam<- read.table("expanded_genes_in_ea_fam.txt", header = FALSE, fill = TRUE, sep = "\t", row.names=1)


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

row.names(expanded_genes_in_ea_fam_formatted)<- rownames(expanded_genes_in_ea_fam)
row.names(contracted_genes_in_ea_fam_formatted)<- rownames(contracted_genes_in_ea_fam)

###print and finish in bash 
#first replace spaces for printing in the GO file's go terms 
GO_files$goName<- GO_files[, gsub(" ", "_", GO_files$goName)]

write.table(expanded_genes_in_ea_fam_formatted, "genes_in_expanded_fams.tab", sep="\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
write.table(contracted_genes_in_ea_fam_formatted, "genes_in_contracted_fams.tab", sep="\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
write.table(GO_files, "GO_files.tab", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
