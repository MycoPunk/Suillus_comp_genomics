#load libraries
library("seqinr")
library("data.table")

options(stringsAsFactors = FALSE)

#read in the hmm results
input<- fread(input = "hmm_output_file", 
              stringsAsFactors=FALSE, 
              header = FALSE, sep = ";")
#remove hashed lines 
input<- input[!grepl("#", input$V1),]

#read in the target
target<- data.frame(t(read.table(file = "one_out.fasta")))
results_df<- target$X1
  
#is the target positive or negative?
 #does it start with >Class or >NEG?
  if(startsWith(target$X1, ">Class")){
    results_df<- cbind(results_df,"GPCR")
  } else {
    results_df<- cbind(results_df,"Non_GPCR")
  }
colnames(results_df)<- c("gene_name", "target_is")
results_df<- data.frame(results_df)

#log the call
#did it match is to a GPCR in the HMM model?
if(nrow(input) <1) {
  results_df<-cbind(results_df,"Non-GPCR")  
} else {
  results_df<-cbind(results_df,"GPCR")
}
colnames(results_df)<- c("gene_name", "target_is", "call_is")
results_df<- data.frame(results_df)

#log false positives
# false positive would be when it is NOT a GPCR, but it gets a positive hit (nrow(input is == 1))'
  if(results_df$target_is =="GPCR") {
    results_df<- cbind(results_df,"-")
  } else if (nrow(input) ==1) {
    results_df<-cbind(results_df,"Y")  
  } else {
    results_df<-cbind(results_df,"N")
}
colnames(results_df)<- c("gene_name", "target_is", "call_is", "false_pos")
results_df<- data.frame(results_df)

#log false negatives
#a false negative is when it IS a GPCR, but it didn't find a match (nrow(input is not == 1))'
if(results_df$target_is =="Non_GPCR") {
  results_df<- cbind(results_df,"-")
} else if (nrow(input) <1) {
  results_df<-cbind(results_df,"Y")  
} else {
  results_df<-cbind(results_df,"N")
}
colnames(results_df)<- c("gene_name", "target_is", "call_is", "false_pos", "false_neg")
results_df<- data.frame(results_df)

#log of take home result (is it a correct call?) 
if((results_df$false_pos != "Y") & (results_df$false_neg != "Y")){
  results_df<- cbind(results_df,"Y")
} else {
  results_df<- cbind(results_df,"N")
}
colnames(results_df)<- c("gene_name","target_is","call_is","false_pos","false_neg","correct_call")
results_df<- data.frame(results_df)

#write it to an output file
write.table(results_df, row.names = FALSE, col.names = FALSE, quote = FALSE, file = "GPCR_itteration_totals.csv", append = TRUE)
