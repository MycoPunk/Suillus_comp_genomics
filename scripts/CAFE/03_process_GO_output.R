setwd("~/bigdata/Suillus_comp_genomics/All_GOterms")
options(stringsAsFactors = FALSE)

#read in files
expansions_files<- list.files(path ="~/bigdata/Suillus_comp_genomics/All_GOterms/", pattern = "GO_terms_output.test.expansions*", full.names=T)
contractions_files<- list.files(path ="~/bigdata/Suillus_comp_genomics/All_GOterms/", pattern = "GO_terms_output.test.contractions*", full.names=T)

#exclude empty files
expansions_files <- expansions_files[which(file.info(expansions_files)$size>0)]
contractions_files <- contractions_files[which(file.info(contractions_files)$size>0)]

expansions<-lapply(expansions_files, fread, sep = "\t", header = FALSE, colClasses = 'character', stringsAsFactors=FALSE)
contractions<-lapply(contractions_files, fread, sep = "\t", header = FALSE, colClasses = 'character', stringsAsFactors=FALSE)

#keep names consistant with gene family #s
expansions_names<- gsub(".*test.expansions_", "", expansions_files)
expansions_names<- gsub(".txt", "", expansions_names)
names(expansions)<- expansions_names

contractions_names<- gsub(".*test.contractions_", "", contractions_files)
contractions_names<- gsub(".txt", "", contractions_names)
names(contractions)<- contractions_names



#expansions<- rbindlist(sapply(expansions_names, fread, simplify = FALSE),
#                      use.names = TRUE, idcol = "FileName", fill = TRUE, header = TRUE)
                      
#contractions<- rbindlist(sapply(contractions_names, fread, simplify = FALSE),
#                       use.names = TRUE, idcol = "FileName", fill = TRUE)


#call table on each data frame 
expansions_table<- lapply(expansions, table)
contractions_table<- lapply(contractions, table)

#these are the gene families of interest (expanded in more than four species)
expansions_table$`289`
expansions_table$`279`
expansions_table$`130`
expansions_table$`197`
expansions_table$`233`
expansions_table$`311`
expansions_table$`168`
expansions_table$`320`
expansions_table$`198`
expansions_table$`233`
expansions_table$`310`



contractions_table$`61`
contractions_table$`108`
contractions_table$`29`
contractions_table$`104`
contractions_table$`185`
contractions_table$`119`
contractions_table$`122`
contractions_table$`248`
contractions_table$`172`
contractions_table$`184`
contractions_table$`197`
contractions_table$`203`
contractions_table$`21`
contractions_table$`66`
contractions_table$`126`
contractions_table$`198`
contractions_table$`247`

contractions_table$`60`
contractions_table$`127`
contractions_table$`147`
contractions_table$`49`
contractions_table$`92`
