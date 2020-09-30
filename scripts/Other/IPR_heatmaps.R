setwd("~/bigdata/Suillus_comp_genomics/All_InterPro")
options(stringsAsFactors = FALSE)

#laod libraries
library('pheatmap')
library('tidyr')
library('dplyr')
library('topGO')

###Interpro2Go mapping

#read in mapping file
Interpro2GO_mapping_file<- read.table("Interpro2GO_mapping_file.txt", sep= "\t", quote = "")
#format to split collumns 
Interpro2GO_mapping_file2<- data.frame(sub(' ', ';', Interpro2GO_mapping_file$V1))
Interpro2GO_mapping_file3<- data.frame(sub('>', ';', Interpro2GO_mapping_file2$sub...........Interpro2GO_mapping_file.V1.))
#split based on ";"
colnames(Interpro2GO_mapping_file3)<- "merged"
Interpro2GO_mapping_file4 <- data.frame(do.call('rbind', strsplit(as.character(Interpro2GO_mapping_file3$merged),';',fixed=TRUE)))
colnames(Interpro2GO_mapping_file4)<- c("InterPro_ID", "InterPro_des", "GO_des", "GO_ID")

#clean 
Interpro2GO_mapping_file4$InterPro_ID_clean<- sub('InterPro:', '', Interpro2GO_mapping_file4$InterPro_ID)

#clean up
rm(Interpro2GO_mapping_file)
rm(Interpro2GO_mapping_file2)
rm(Interpro2GO_mapping_file3)

#read in interpro files
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

#call table to get totals
Interpro_table<- data.frame(table(Interpro_files3$FileName, Interpro_files3$iprId), stringsAsFactors = FALSE)
names(Interpro_table)<- c("Species", "iprId", "value")

#clean up
rm(Interpro_files)
rm(Interpro_files2)
rm(Interpro_files3)

input<- Interpro2GO_mapping_file4 %>% 
  group_by(InterPro_ID_clean)%>%
  summarise(GO_ID = paste(unique(GO_ID), collapse=','))

#print and read it in
write.table(input, "IPR_to_go_mapping_formatted.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
geneID2GO <- readMappings(file = "IPR_to_go_mapping_formatted.txt") 

##set gene universe
geneUniverse <- names(geneID2GO) 

##get genes of interest 
#put into wide form (this takes a hot minute)
Interpro_table_reshaped<- reshape(Interpro_table, idvar = "Species", timevar = "iprId", direction = "wide")
#fix the names
names(Interpro_table_reshaped)<- sub("value.", "", names(Interpro_table_reshaped))

#turn first col. into row names and remove the first col. 
row.names(Interpro_table_reshaped)<- Interpro_table_reshaped[,1]
Interpro_table_reshaped<- Interpro_table_reshaped[,-1]
#split the datasets
Interpro_for_rendering_Suillus<- Interpro_table_reshaped[grep("Sui", row.names(Interpro_table_reshaped)),]
Interpro_for_rendering_Other<- Interpro_table_reshaped[grep("Sui", row.names(Interpro_table_reshaped), invert = TRUE),]
#remove suspects
suspects<- c("Gyrli1", "Rhivul1", "Rhivi1", "Rhives1", "Rhisa1")
Interpro_for_rendering_Other_wo<- Interpro_for_rendering_Other[!rownames(Interpro_for_rendering_Other) %in% suspects, ]

##get significantly different domains
#Rewrite to get only the domains significantly GREATER in Suillus
my.t.test.p.value <- function(...) {
  obj<-try(t.test(..., alternative = 'greater'), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

#run the function in place of the t-test call. 
t_test_results_greater<- Map(my.t.test.p.value,Interpro_for_rendering_Suillus,Interpro_for_rendering_Other)
length(t_test_results_greater[t_test_results_greater < .001])

#remove NAs
length(t_test_results_greater[t_test_results_greater == "NA"])
t_test_results_greater_clean<- t_test_results_greater[t_test_results_greater != "NA"]
length(t_test_results_greater)
length(t_test_results_greater_clean)
length(t_test_results_greater_clean[t_test_results_greater_clean == "NA"])
#subset list
t_test_results_greater_list<- t_test_results_greater_clean[t_test_results_greater_clean < .001]
length(t_test_results_greater_list)

#with High and Moderate removed
#run the function in place of the t-test call. 
t_test_results_greater_wo<- Map(my.t.test.p.value,Interpro_for_rendering_Suillus,Interpro_for_rendering_Other_wo)
length(t_test_results_greater_wo[t_test_results_greater_wo < .001])

#remove NAs
length(t_test_results_greater_wo[t_test_results_greater_wo == "NA"])
t_test_results_greater_wo_clean<- t_test_results_greater_wo[t_test_results_greater_wo != "NA"]
length(t_test_results_greater_wo)
length(t_test_results_greater_wo_clean)
length(t_test_results_greater_wo_clean[t_test_results_greater_wo_clean == "NA"])
#subset list
t_test_results_greater_wo_list<- t_test_results_greater_wo_clean[t_test_results_greater_wo_clean < .001]
length(t_test_results_greater_wo_list)


#Rewrite to get only the domains significantly LESS in Suillus
my.t.test.p.value <- function(...) {
  obj<-try(t.test(..., alternative = 'less'), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

#run the function in place of the t-test call. 
t_test_results_less<- Map(my.t.test.p.value,Interpro_for_rendering_Suillus,Interpro_for_rendering_Other)
length(t_test_results_less[t_test_results_less < .001])
#remove NAs
length(t_test_results_less[t_test_results_less == "NA"])
t_test_results_less_clean<- t_test_results_less[t_test_results_less != "NA"]
length(t_test_results_less)
length(t_test_results_less_clean)
length(t_test_results_less_clean[t_test_results_less_clean == "NA"])

#subset list
t_test_results_less_list<- t_test_results_less_clean[t_test_results_less_clean < .001]
length(t_test_results_less_list)


#for High and Moderate removed
#run the function in place of the t-test call. 
t_test_results_less_wo<- Map(my.t.test.p.value,Interpro_for_rendering_Suillus,Interpro_for_rendering_Other_wo)
length(t_test_results_less_wo[t_test_results_less_wo < .001])
#remove NAs
length(t_test_results_less_wo[t_test_results_less_wo == "NA"])
t_test_results_less_wo_clean<- t_test_results_less_wo[t_test_results_less_wo != "NA"]
length(t_test_results_less_wo)
length(t_test_results_less_wo_clean)
length(t_test_results_less_wo_clean[t_test_results_less_wo_clean == "NA"])

#subset list
t_test_results_less_wo_list<- t_test_results_less_wo_clean[t_test_results_less_wo_clean < .001]
length(t_test_results_less_wo_list)



#set genes of interest 
genesOfInterest_greater<- names(t_test_results_greater_list)
genesOfInterest_less<- names(t_test_results_less_list)
genesOfInterest_greater_wo<- names(t_test_results_greater_wo_list)
genesOfInterest_less_wo<- names(t_test_results_less_wo_list)

#tell topGO where to look for the genes of interest
geneList_greater <- factor(as.integer(geneUniverse %in% genesOfInterest_greater))
geneList_less <- factor(as.integer(geneUniverse %in% genesOfInterest_less))
geneList_greater_wo <- factor(as.integer(geneUniverse %in% genesOfInterest_greater_wo))
geneList_less_wo <- factor(as.integer(geneUniverse %in% genesOfInterest_less_wo))
names(geneList_greater) <- geneUniverse
names(geneList_less) <- geneUniverse
names(geneList_greater_wo) <- geneUniverse
names(geneList_less_wo) <- geneUniverse


#to put the data into an object of type 'topGOdata'. 
#This will contain the list of genes of interest, the GO annotations, and the GO hierarchy.
myGOdata_greater <- new("topGOdata", description="greater", ontology="BP", allGenes=geneList_greater,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
myGOdata_less <- new("topGOdata", description="less", ontology="BP", allGenes=geneList_less,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
myGOdata_greater_wo <- new("topGOdata", description="greater", ontology="BP", allGenes=geneList_greater_wo,  annot = annFUN.gene2GO, gene2GO = geneID2GO)
myGOdata_less_wo <- new("topGOdata", description="less", ontology="BP", allGenes=geneList_less_wo,  annot = annFUN.gene2GO, gene2GO = geneID2GO)


#NOTE The 'ontology' argument can be set to 'BP' (biological process), 'MF' (molecular function), or 'CC' (cellular component).
#NOTE: significant genes are just what topGO calls your genes of interest that had GO mapping
myGOdata_greater
myGOdata_less
myGOdata_greater_wo
myGOdata_less_wo

#NOTE access genes of interest like this
#sg <- sigGenes(myGOdata_greater)
#str(sg)
#numSigGenes(myGOdata_greater)

#test for enrichment 
resultFisher_greater <- runTest(myGOdata_greater, algorithm="weight01", statistic="fisher")
resultFisher_greater

resultFisher_less <- runTest(myGOdata_less, algorithm="weight01", statistic="fisher")
resultFisher_less

resultFisher_greater_wo <- runTest(myGOdata_greater_wo, algorithm="weight01", statistic="fisher")
resultFisher_greater_wo

resultFisher_less_wo <- runTest(myGOdata_less_wo, algorithm="weight01", statistic="fisher")
resultFisher_less_wo

#get top 10 results
allRes_greater <- GenTable(myGOdata_greater, classicFisher = resultFisher_greater, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 20)
allRes_greater
allRes_greater_df<- as.data.frame(allRes_greater)
#18 results significant at p<0.01
#GO.ID                                        Term Annotated Significant Expected classicFisher
#1  GO:0055114                 oxidation-reduction process       924          65    26.21       7.1e-13
#2  GO:0006556 S-adenosylmethionine biosynthetic proces...         6           5     0.17       1.0e-07
#3  GO:0006412                                 translation       444          36    12.60       4.2e-07
#4  GO:0009113      purine nucleobase biosynthetic process         7           5     0.20       9.0e-06
#5  GO:0006419                  alanyl-tRNA aminoacylation         4           3     0.11       8.9e-05
#6  GO:0009058                        biosynthetic process      2845         109    80.71       0.00042
#7  GO:0000256                 allantoin catabolic process         6           3     0.17       0.00042
#8  GO:0000154                           rRNA modification        14           3     0.40       0.00080
#9  GO:0006481              C-terminal protein methylation         2           2     0.06       0.00080
#10 GO:0006099                    tricarboxylic acid cycle        27           5     0.77       0.00086
#11 GO:0055085                     transmembrane transport       422          23    11.97       0.00124
#12 GO:0005975              carbohydrate metabolic process       490          21    13.90       0.00133
#13 GO:0006006                   glucose metabolic process        32           4     0.91       0.00233
#14 GO:0019551 glutamate catabolic process to 2-oxoglut...         3           2     0.09       0.00236
#15 GO:0001522                     pseudouridine synthesis        21           4     0.60       0.00260
#16 GO:0019915                               lipid storage         4           2     0.11       0.00464
#17 GO:0006097                            glyoxylate cycle         5           2     0.14       0.00758
#18 GO:0006891      intra-Golgi vesicle-mediated transport         5           2     0.14       0.00758

allRes_less <- GenTable(myGOdata_less, classicFisher = resultFisher_less, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
allRes_less
allRes_less_df<- as.data.frame(allRes_less)
#three significant esults at p <0.01
#GO.ID                                        Term Annotated Significant Expected classicFisher
#  GO:0006412                                 translation       444          18     4.68       1.4e-06
#  GO:0000103                        sulfate assimilation         8           2     0.08        0.0030
#  GO:0045454                      cell redox homeostasis        10           2     0.11        0.0047


allRes_greater_wo <- GenTable(myGOdata_greater_wo, classicFisher = resultFisher_greater_wo, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 20)
allRes_greater_wo
allRes_greater_wo_df<- as.data.frame(allRes_greater_wo)
#there are 9 significant results
#        GO.ID                                        Term Annotated Significant Expected classicFisher
#1  GO:0055114                 oxidation-reduction process       924          55    18.89       6.1e-14
#2  GO:0009058                        biosynthetic process      2845          70    58.17       1.0e-04
#3  GO:0000256                 allantoin catabolic process         6           3     0.12       0.00016
#4  GO:0000154                           rRNA modification        14           3     0.29       0.00041
#5  GO:0055085                     transmembrane transport       422          19     8.63       0.00078
#6  GO:0019551 glutamate catabolic process to 2-oxoglut...         3           2     0.06       0.00123
#7  GO:0006099                    tricarboxylic acid cycle        27           4     0.55       0.00206
#8  GO:0006097                            glyoxylate cycle         5           2     0.10       0.00400
#9  GO:0001522                     pseudouridine synthesis        21           3     0.43       0.00855

allRes_less_wo <- GenTable(myGOdata_less_wo, classicFisher = resultFisher_less_wo, orderBy = "resultFisher", ranksOf = "classicFisher", topNodes = 10)
allRes_less_wo
allRes_less_wo_df<- as.data.frame(allRes_less_wo)
#GO.ID                                        Term Annotated Significant Expected classicFisher
#1  GO:0045454                      cell redox homeostasis        10           2     0.05        0.0012


##link the GO terms back to the IPR annotations
#read back in files
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

##get IPR Ids associated with GO terms #GO:0055114 = oxidation-reduction process (greater) 
#subset to only the two relevant GO terms 
Interpro2GO_mapping_file4_greater<- Interpro2GO_mapping_file4[Interpro2GO_mapping_file4$GO_ID == " GO:0055114",]

#shrink the input dataframe to only the relevent IPR terms under those GO terms
shrunk_input_greater<- Interpro_files3[Interpro_files3$iprId %in% Interpro2GO_mapping_file4_greater$InterPro_ID_clean,]
nrow(Interpro_files3)
nrow(shrunk_input_greater)

##now use the above as your input 
#call table to get totals
Interpro_table_greater<- data.frame(table(shrunk_input_greater$FileName, shrunk_input_greater$iprDesc), stringsAsFactors = FALSE)
names(Interpro_table_greater)<- c("Species", "iprDesc", "value")

#put into wide form (this takes a hot minute)
Interpro_table_reshaped_greater<- reshape(Interpro_table_greater, idvar = "Species", timevar = "iprDesc", direction = "wide")

#fix the names
names(Interpro_table_reshaped_greater)<- sub("value.", "", names(Interpro_table_reshaped_greater))

#turn first col. into row names and remove the first col. 
row.names(Interpro_table_reshaped_greater)<- Interpro_table_reshaped_greater[,1]
Interpro_table_reshaped_greater<- Interpro_table_reshaped_greater[,-1]

#clean up 
rm(Interpro_files)
rm(Interpro_files2)
rm(Interpro_files3)
rm(Interpro_table)

#before rendering - Move Thega1 and Theter1, so that they don't end up below suillus
#names in order
names_in_order<- c("Suibov1",
                   "Suibr1",
                   "Suicot1",
                   "Suidec1",
                   "Suifus1",
                   "Suihi1",
                   "Suilu4",
                   "Suiocc1",
                   "Suisu1",
                   "Suitom1",
                   "Suivar1",
                   "Suiame1",
                   "Suidis1",
                   "Suipic1",
                   "Suipla1",
                   "Suiplo1",
                   "Suisubl1",
                   "Suiamp1",
                   "Suicli1",
                   "Suipal1",
                   "Suisub1",
                   "Suigr1",
                   "Suilak1",
                   "Rhisa1",
                   "Rhitru1",
                   "Rhives1",
                   "Rhivi1",
                   "Rhivul1",
                   "Amamu1",
                   "Cananz1",
                   "Gaumor1_1",
                   "Gyrli1",
                   "Hebcy2",
                   "Hydru2",
                   "Hyssto1",
                   "Lacam2",
                   "Lacbi2",
                   "Paxin1",
                   "Pilcr1",
                   "Pismi1",
                   "Pisti1",
                   "Rusbre1",
                   "Ruscom1",
                   "Sclci1",
                   "Theter1",
                   "Thega1")


Interpro_table_reshaped_by_greater<- Interpro_table_reshaped_greater[match(names_in_order, rownames(Interpro_table_reshaped_greater)),]                 

#split the datasets
Interpro_for_rendering_Suillus_greater<- Interpro_table_reshaped_by_greater[grep("Sui", row.names(Interpro_table_reshaped_by_greater)),]
Interpro_for_rendering_Other_greater<- Interpro_table_reshaped_by_greater[grep("Sui", row.names(Interpro_table_reshaped_by_greater), invert = TRUE),]

#function to get only the domains significantly GREATER in Suillus
my.t.test.p.value <- function(...) {
  obj<-try(t.test(..., alternative = 'greater'), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

#run the function in place of the t-test call. 
t_test_results_greater<- Map(my.t.test.p.value,Interpro_for_rendering_Suillus_greater,Interpro_for_rendering_Other_greater)
length(t_test_results_greater[t_test_results_greater < .001])

length(t_test_results_greater[t_test_results_greater == "NA"])
t_test_results_greater_clean<- t_test_results_greater[t_test_results_greater != "NA"]
length(t_test_results_greater)
length(t_test_results_greater_clean)
length(t_test_results_greater_clean[t_test_results_greater_clean == "NA"])

#subset list
t_test_results_greater_list<- t_test_results_greater_clean[t_test_results_greater_clean < .001]
length(t_test_results_greater_list)
t_test_results_greater_list_in_total<- t_test_results_greater_clean[t_test_results_greater_clean < .001]
length(t_test_results_greater_list_in_total)

#put p.vals in order 
t_test_results_greater_list_to_sort<-unlist(t_test_results_greater_list)
t_test_results_greater_list_df<- as.data.frame(t_test_results_greater_list_to_sort)
t_test_results_greater_list_ordered<- t_test_results_greater_list_df[order(t_test_results_greater_list_df$t_test_results_greater_list_to_sort), , drop = FALSE]

#subset dataframes for plotting
Interpro_for_rendering_greater<- Interpro_table_reshaped_by_greater[,names(Interpro_table_reshaped_by_greater) %in% rownames(t_test_results_greater_list_ordered)]
ncol(Interpro_for_rendering_greater)

Interpro_for_rendering_greater_all<- Interpro_table_reshaped_by_greater[,names(Interpro_table_reshaped_by_greater) %in% names(t_test_results_greater_list_in_total)]
ncol(Interpro_for_rendering_greater_all)


#print for supplemental table 
write.table(Interpro_for_rendering_greater_all, "IPR_redox.txt", sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
dim(Interpro_for_rendering_greater_all)
###plot this: 

##make color gradient 
library(RColorBrewer)
library(viridis)
pheatmap(Interpro_for_rendering_greater[,1:ncol(Interpro_for_rendering_greater)], 
         angle_col = "315", 
         fontsize = 3,
         cluster_rows= FALSE,
         cluster_cols = TRUE,
         scale = "column",
         color = inferno(100),
         #color = colfunc(100),
         #breaks = 1,
         border_color = NA,
         cellheight=5,cellwidth=5,
         #display_numbers = Interpro_for_rendering_greater,
         filename = "heatmap_IPR_ox_redux.pdf")


###look at highest differentiated p-fam domains overall (not structured by highest ranked GO annotation)
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

#call table to get totals
Interpro_table<- data.frame(table(Interpro_files3$FileName, Interpro_files3$iprDesc), stringsAsFactors = FALSE)
names(Interpro_table)<- c("Species", "iprDesc", "value")

#put into wide form (this takes a hot minute)
Interpro_table_reshaped<- reshape(Interpro_table, idvar = "Species", timevar = "iprDesc", direction = "wide")

#fix the names
names(Interpro_table_reshaped)<- sub("value.", "", names(Interpro_table_reshaped))

#turn first col. into row names and remove the first col. 
row.names(Interpro_table_reshaped)<- Interpro_table_reshaped[,1]
Interpro_table_reshaped<- Interpro_table_reshaped[,-1]

#reorder using the names ist defined in section 1
Interpro_table_reshaped_by_name<- Interpro_table_reshaped[match(names_in_order, rownames(Interpro_table_reshaped)),]                 
Interpro_table_reshaped<- Interpro_table_reshaped_by_name

#split the datasets
Interpro_for_rendering_Suillus<- Interpro_table_reshaped[grep("Sui", row.names(Interpro_table_reshaped)),]
Interpro_for_rendering_Other<- Interpro_table_reshaped[grep("Sui", row.names(Interpro_table_reshaped), invert = TRUE),]

##get significantly different domains
#Rewrite to get only the domains significantly GREATER in Suillus
my.t.test.p.value <- function(...) {
  obj<-try(t.test(..., alternative = 'greater'), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

#run the function in place of the t-test call. 
t_test_results_greater<- Map(my.t.test.p.value,Interpro_for_rendering_Suillus,Interpro_for_rendering_Other)
length(t_test_results_greater[t_test_results_greater < .001])

#remove NAs
length(t_test_results_greater[t_test_results_greater == "NA"])
t_test_results_greater_clean<- t_test_results_greater[t_test_results_greater != "NA"]
length(t_test_results_greater)
length(t_test_results_greater_clean)
length(t_test_results_greater_clean[t_test_results_greater_clean == "NA"])
#subset list
t_test_results_greater_list<- t_test_results_greater_clean[t_test_results_greater_clean < .001]
length(t_test_results_greater_list)


#get only the domains significantly LESS in Suillus
my.t.test.p.value <- function(...) {
  obj<-try(t.test(..., alternative = 'less'), silent=TRUE)
  if (is(obj, "try-error")) return(NA) else return(obj$p.value)
}

#run the function in place of the t-test call. 
t_test_results_less<- Map(my.t.test.p.value,Interpro_for_rendering_Suillus,Interpro_for_rendering_Other)
length(t_test_results_less[t_test_results_less < .001])
#remove NAs
length(t_test_results_less[t_test_results_less == "NA"])
t_test_results_less_clean<- t_test_results_less[t_test_results_less != "NA"]
length(t_test_results_less)
length(t_test_results_less_clean)
length(t_test_results_less_clean[t_test_results_less_clean == "NA"])

#subset list
t_test_results_less_list<- t_test_results_less_clean[t_test_results_less_clean < .001]
length(t_test_results_less_list)

#put p.vals in order 
t_test_results_greater_list_to_sort<-unlist(t_test_results_greater_list)
t_test_results_greater_list_df<- as.data.frame(t_test_results_greater_list_to_sort)
t_test_results_greater_list_ordered<- t_test_results_greater_list_df[order(t_test_results_greater_list_df$t_test_results_greater_list_to_sort), , drop = FALSE]

#put p.vals in order 
t_test_results_less_list_to_sort<-unlist(t_test_results_less_list)
t_test_results_less_list_df<- as.data.frame(t_test_results_less_list_to_sort)
t_test_results_less_list_ordered<- t_test_results_less_list_df[order(t_test_results_less_list_df$t_test_results_less_list_to_sort), , drop = FALSE]

#get top 90 values 
t_test_results_greater_list_ordered_top50<- t_test_results_greater_list_ordered[1:90,, drop = FALSE]
t_test_results_less_list_ordered_top50<- t_test_results_less_list_ordered[1:90,, drop = FALSE]

#subset dataframes for plotting
Interpro_for_rendering_greater<- Interpro_table_reshaped[,names(Interpro_table_reshaped) %in% rownames(t_test_results_greater_list_ordered_top50)]
ncol(Interpro_for_rendering_greater)

Interpro_for_rendering_less<- Interpro_table_reshaped[,names(Interpro_table_reshaped) %in% rownames(t_test_results_less_list_ordered_top50)]
ncol(Interpro_for_rendering_less)

#subset dataframes for Supplemental tables 
Interpro_for_rendering_greater_suptable<- Interpro_table_reshaped[,names(Interpro_table_reshaped) %in% rownames(t_test_results_greater_list_ordered)]
nrow(t_test_results_greater_list_ordered)
ncol(Interpro_for_rendering_greater_suptable)

Interpro_for_rendering_less_suptable<- Interpro_table_reshaped[,names(Interpro_table_reshaped) %in% rownames(t_test_results_less_list_ordered)]
nrow(t_test_results_less_list_ordered)
ncol(Interpro_for_rendering_less_suptable)

##plot both of these: 
###make color gradient 
pheatmap(Interpro_for_rendering_greater[,1:ncol(Interpro_for_rendering_greater)], 
         angle_col = "315", 
         fontsize = 3,
         cluster_rows= FALSE,
         cluster_cols = TRUE,
         scale = "column",
         color = inferno(100),
         #color = colfunc(100),
         #breaks = 1,
         border_color = NA,
         cellheight=5,cellwidth=5,
         #display_numbers = Interpro_for_rendering_greater,
         filename = "heatmap_PIR_square_greater.pdf")


#write file for Sup table
write.table(Interpro_for_rendering_greater_suptable, "IPR_over_rep_in_Suillus.txt", sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)


pheatmap(Interpro_for_rendering_less[,1:ncol(Interpro_for_rendering_less)], 
         angle_col = "315", 
         fontsize = 3,
         cluster_rows= FALSE,
         cluster_cols = TRUE,
         scale = "column",
         color = inferno(100),
         #color = colfunc(100),
         #breaks = 1,
         border_color = NA,
         cellheight=5,cellwidth=5,
         #display_numbers = Interpro_for_rendering_less,
         filename = "heatmap_PIR_square_less.pdf")

#write file for Sup table
write.table(Interpro_for_rendering_less_suptable, "IPR_under_rep_in_Suillus.txt", sep="\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
