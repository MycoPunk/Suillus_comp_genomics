#load libraries
library("seqinr")
library("data.table")
library("stringr")
library("tidyr")
library("dplyr")

options(stringsAsFactors = FALSE)
setwd("~/bigdata/GPCR/HMM_model_seqs/expanded_GPCR_set/expanded_fastas/Class1")
#read in the input file
expanded_set<- seqinr::read.fasta(file = "class1_expanded.fasta", 
                                  seqtype = "AA",as.string = TRUE, set.attributes = FALSE)

Positive_set<- seqinr::read.fasta(file = "GPCR_class1_whole.fasta", 
                                  seqtype = "AA",as.string = TRUE, set.attributes = FALSE)

Negative_set<- seqinr::read.fasta(file = "non-GPCR_whole.fasta", 
                                  seqtype = "AA",as.string = TRUE, set.attributes = FALSE)

#format all input files
Positive_set_df<- t(data.frame(Positive_set))
colnames(Positive_set_df)<- "seqs"
rownames(Positive_set_df) <- paste0(">",  rownames(Positive_set_df))

expanded_set_df<- t(data.frame(expanded_set))
colnames(expanded_set_df)<- "seqs"
rownames(expanded_set_df) <- paste0(">",  rownames(expanded_set_df))

#subset the negative input file to the same number to files as the positive inputs 
n_neg<- nrow(Positive_set_df)
negative_set_equal_size <- Negative_set[ sample(seq_along(Negative_set),size=n_neg,replace=FALSE)]
#and format the negative set like the rest of the sets

Negative_set_df<- t(data.frame(negative_set_equal_size))
colnames(Negative_set_df)<- "seqs"
rownames(Negative_set_df) <- paste0(">",  rownames(Negative_set_df))

#combine the positive and negative sets into one set to draw from
set_todraw_from<- data.frame(rbind(Negative_set_df, Positive_set_df))
set_todraw_from_df<- data.frame(cbind(row.names(set_todraw_from), set_todraw_from))
colnames(set_todraw_from_df)<- c("names", "seqs")

#randomly sample one fasta from the set
one_fasta <- data.frame(sample_n(set_todraw_from_df, 1, row.names = TRUE))

#remove the one_fasta from the list of positives (in case it's in there)
minus_one_fasta<- data.frame(Positive_set_df[!(row.names(Positive_set_df) %in% one_fasta$names), ])
colnames(minus_one_fasta)<- c("seqs")

#add the remaining positive dataset to the expanded dataset (found by psiblasting the positive dataset)
set_for_model<- rbind(expanded_set_df, minus_one_fasta)

minus_one_fasta<- data.frame(cbind(row.names(set_for_model), set_for_model))
colnames(minus_one_fasta)<- c("names", "seqs")

#interleaf the seqs and ID's
format<- do.call(rbind, lapply(seq(nrow(minus_one_fasta)), function(i) t(minus_one_fasta[i, ])))
format2<- do.call(rbind, lapply(seq(nrow(one_fasta)), function(i) t(one_fasta[i, ])))


#print the new test fasta and the -1 dataset to build in clustal and test the HMM 
write.table(format, row.names = FALSE, col.names = FALSE, quote = FALSE, file = "minus_one.fasta")
write.table(format2, row.names = FALSE, col.names = FALSE, quote = FALSE, file = "one_out.fasta")
