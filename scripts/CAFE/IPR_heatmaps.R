setwd("~/bigdata/Suillus_comp_genomics/All_InterPro")
options(stringsAsFactors = FALSE)

#laod libraries
library('pheatmap')

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

#fix empties
Interpro_table_reshaped_filtered <- Interpro_table_reshaped +.001

#clean up 
rm(Interpro_files)
rm(Interpro_files2)
rm(Interpro_files3)
rm(Interpro_table)
rm(Interpro_table_reshaped)

ncol(Interpro_table_reshaped_filtered)
#there are 7960 gene fams - you'e going to have to split this up somehow. 

#before rendering - Move Thega1 and Theter1, so that they don't end up below suillus
nrow(Interpro_table_reshaped_filtered)
#short<- Interpro_table_reshaped_filtered[1:44,]
last.two<- Interpro_table_reshaped_filtered[45:46,]
others<- Interpro_table_reshaped_filtered[1:21,]
Suilluses<- Interpro_table_reshaped_filtered[22:44,]
Interpro_for_rendering<- rbind(Suilluses, others, last.two)

#print all interpro heatmaps in sections for visual inspection
pdf("pheatmap1.pdf", width=12, height=9)
out <- pheatmap(Interpro_for_rendering[,1:200], 
                show_rownames=T,
                #cluster_cols=T,
                cluster_rows=F,
                #scale="row",
                border_color=FALSE,
                cex=0.7)
out2 <- pheatmap(Interpro_for_rendering[,200:400], 
                show_rownames=T,
                #cluster_cols=T,
                cluster_rows=F,
                #scale="row",
                border_color=FALSE,
                cex=0.7)
out3 <- pheatmap(Interpro_for_rendering[,400:600], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out4 <- pheatmap(Interpro_for_rendering[,600:800], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out5 <- pheatmap(Interpro_for_rendering[,800:1000], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
dev.off()


pdf("pheatmap2.pdf", width=12, height=9)
out <- pheatmap(Interpro_for_rendering[,1000:1200], 
                show_rownames=T,
                #cluster_cols=T,
                cluster_rows=F,
                #scale="row",
                border_color=FALSE,
                cex=0.7)
out2 <- pheatmap(Interpro_for_rendering[,1200:1400], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out3 <- pheatmap(Interpro_for_rendering[,1400:1600], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out4 <- pheatmap(Interpro_for_rendering[,1600:1800], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out5 <- pheatmap(Interpro_for_rendering[,1800:2000], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
dev.off()

pdf("pheatmap3.pdf", width=12, height=9)
out <- pheatmap(Interpro_for_rendering[,2000:2200], 
                show_rownames=T,
                #cluster_cols=T,
                cluster_rows=F,
                #scale="row",
                border_color=FALSE,
                cex=0.7)
out2 <- pheatmap(Interpro_for_rendering[,2200:2400], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out3 <- pheatmap(Interpro_for_rendering[,2400:2600], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out4 <- pheatmap(Interpro_for_rendering[,2600:2800], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out5 <- pheatmap(Interpro_for_rendering[,2800:3000], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
dev.off()

pdf("pheatmap4.pdf", width=12, height=9)
out <- pheatmap(Interpro_for_rendering[,3000:3200], 
                show_rownames=T,
                #cluster_cols=T,
                cluster_rows=F,
                #scale="row",
                border_color=FALSE,
                cex=0.7)
out2 <- pheatmap(Interpro_for_rendering[,3200:3400], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out3 <- pheatmap(Interpro_for_rendering[,3400:3600], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out4 <- pheatmap(Interpro_for_rendering[,3600:3800], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out5 <- pheatmap(Interpro_for_rendering[,3800:4000], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
dev.off()

pdf("pheatmap5.pdf", width=12, height=9)
out <- pheatmap(Interpro_for_rendering[,4000:4200], 
                show_rownames=T,
                #cluster_cols=T,
                cluster_rows=F,
                #scale="row",
                border_color=FALSE,
                cex=0.7)
out2 <- pheatmap(Interpro_for_rendering[,4200:4400], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out3 <- pheatmap(Interpro_for_rendering[,4400:4600], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out4 <- pheatmap(Interpro_for_rendering[,4600:4800], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out5 <- pheatmap(Interpro_for_rendering[,4800:5000], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
dev.off()

pdf("pheatmap6.pdf", width=12, height=9)
out <- pheatmap(Interpro_for_rendering[,5000:5200], 
                show_rownames=T,
                #cluster_cols=T,
                cluster_rows=F,
                #scale="row",
                border_color=FALSE,
                cex=0.7)
out2 <- pheatmap(Interpro_for_rendering[,5200:5400], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out3 <- pheatmap(Interpro_for_rendering[,5400:5600], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out4 <- pheatmap(Interpro_for_rendering[,5600:5800], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out5 <- pheatmap(Interpro_for_rendering[,5800:6000], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
dev.off()

pdf("pheatmap7.pdf", width=12, height=9)
out <- pheatmap(Interpro_for_rendering[,6000:6200], 
                show_rownames=T,
                #cluster_cols=T,
                cluster_rows=F,
                #scale="row",
                border_color=FALSE,
                cex=0.7)
out2 <- pheatmap(Interpro_for_rendering[,6200:6400], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out3 <- pheatmap(Interpro_for_rendering[,6400:6600], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out4 <- pheatmap(Interpro_for_rendering[,6600:6800], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out5 <- pheatmap(Interpro_for_rendering[,6800:7000], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
dev.off()


pdf("pheatmap8.pdf", width=12, height=9)
out <- pheatmap(Interpro_for_rendering[,7000:7200], 
                show_rownames=T,
                #cluster_cols=T,
                cluster_rows=F,
                #scale="row",
                border_color=FALSE,
                cex=0.7)
out2 <- pheatmap(Interpro_for_rendering[,7200:7400], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out3 <- pheatmap(Interpro_for_rendering[,7400:7600], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out4 <- pheatmap(Interpro_for_rendering[,7600:7800], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
out5 <- pheatmap(Interpro_for_rendering[,7800:7969], 
                 show_rownames=T,
                 #cluster_cols=T,
                 cluster_rows=F,
                 #scale="row",
                 border_color=FALSE,
                 cex=0.7)
dev.off()


#get only terpene related domains 
terpene_doms<- Interpro_for_rendering[, grepl("Terpene", names(Interpro_for_rendering))]

pheatmap(terpene_doms, 
         angle_col = "45", 
         fontsize = 6,
         cluster_rows= FALSE,
         filename = "terpene_doms.pdf")



###sort the whole df by variance to get at just the interesting IPRs
# get variance for ea col.

#var_list<- apply(Interpro_for_rendering, 2, var)
#the above is variance across all - lets get teh variance between the two groups

Interpro_for_rendering_Suillus<- Interpro_for_rendering[grep("Sui", row.names(Interpro_for_rendering)),]
Interpro_for_rendering_Other<- Interpro_for_rendering[grep("Sui", row.names(Interpro_for_rendering), invert = TRUE),]

#get variance between the two groups 
var_func <- function(x, y) var(x,y)
var_list<- mapply(var_func, Interpro_for_rendering_Suillus, Interpro_for_rendering_Other)


#shrink list to only var > 25
var_list_sig<- var_list[var_list > 25]

length(var_list)
length(var_list_sig)

#subset the main data frame based on inclusion in the subset list
Interpro_for_rendering_subset<- Interpro_for_rendering[,names(Interpro_for_rendering) %in% names(var_list_sig)]

pheatmap(Interpro_for_rendering_subset[,1:ncol(Interpro_for_rendering_subset)], 
         #angle_col = "45", 
         fontsize = 5,
         cluster_rows= FALSE,
         filename = "NEWtest_var25.pdf")


#can we subset this further by significantly different groups? (do a t-test across all for example?)

#shrink list to only var > .0001 to remove any idential values (t-test will throw an error if any of the comparasons are identical)

#shrink list to only var > 25
var_list_sig<- var_list[var_list > 0]

length(var_list)
length(var_list_sig)

#subset the main data frame based on inclusion in the subset list
Interpro_for_rendering_subset<- Interpro_for_rendering[,names(Interpro_for_rendering) %in% names(var_list_sig)]
ncol(Interpro_for_rendering_subset)

#now split the data frame again 
Interpro_for_rendering_Suillus2<- Interpro_for_rendering_subset[grep("Sui", row.names(Interpro_for_rendering_subset)),]
Interpro_for_rendering_Other2<- Interpro_for_rendering_subset[grep("Sui", row.names(Interpro_for_rendering_subset), invert = TRUE),]

#now run a t-test on the subset data frame, where we know the difference is not zero.
t.test_results <- mapply(t.test, x= Interpro_for_rendering_Suillus2, y = Interpro_for_rendering_Other2, SIMPLIFY = F)
ttest.pval <- sapply(t.test_results, '[[', 'p.value')

#reduce list to only significant p-vals
ttest.pval_sig<- ttest.pval[ttest.pval < .00001] # = 1.0e-5
length(ttest.pval_sig)

#now subset the original data frame based on significance
Interpro_for_rendering_subsetNEW<- Interpro_for_rendering[,names(Interpro_for_rendering) %in% names(ttest.pval_sig)]
ncol(Interpro_for_rendering_subsetNEW)
###make color gradient 
#colfunc <- colorRampPalette(c("#432B46", "#FFFFFF"))
library(RColorBrewer)
library(viridis)
pheatmap(Interpro_for_rendering_subsetNEW[,1:ncol(Interpro_for_rendering_subsetNEW)], 
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
         filename = "heatmap_PIR_square_arranged.pdf")

#get names for the significantly different clusters
names_list<- names(Interpro_for_rendering_subsetNEW)
write.table(names_list, "names_of_most_sig_different_interpro_doms.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

help(pheatmap)
