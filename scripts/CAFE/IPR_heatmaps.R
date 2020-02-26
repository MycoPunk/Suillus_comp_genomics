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

#remove factors
#Interpro_table_reshaped_test <- cbind(Interpro_table_reshaped[,1], data.frame(lapply(Interpro_table_reshaped[,2:ncol(Interpro_table_reshaped)], as.numeric), stringsAsFactors=FALSE))

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

#before rendering - rename Thega1 and Theter1, so that they don't end up below suillus
#first try just movig it
nrow(Interpro_table_reshaped_filtered)
short<- Interpro_table_reshaped_filtered[1:44,]
Interpro_for_rendering<- rbind(Interpro_table_reshaped_filtered[45:46,], short)

##try just the first 50
#pheatmap(Interpro_for_rendering[,1:50], 
#         angle_col = "45", 
#         fontsize = 6,
#         cluster_rows= FALSE,
#         filename = "test2.pdf")


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

