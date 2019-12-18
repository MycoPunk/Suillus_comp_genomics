#set wd
setwd("~/bigdata/GPCR/HMM_model_seqs/expanded_GPCR_set/expanded_fastas/seqs/")

#read in libraries
library(dplyr)
library(readr)
library(data.table)
library(seqinr)
library(reshape)
library(ggplot2)

#read in files
filesnames<- list.files(path ="~/bigdata/GPCR/HMM_model_seqs/expanded_GPCR_set/expanded_fastas/seqs/", pattern = "*domtblout*", full.names=T)

HMM_results<-lapply(filesnames, fread, sep = " ", header = FALSE, stringsAsFactors=FALSE)
#HMM_results<-lapply(filesnames, fread, sep = " ", header = FALSE, stringsAsFactors=FALSE)

col_names<- c("target name", 
              "accession", 
              "tlen",
              "query name", 
              "accession", 
              "qlen", 
              "E-value", 
              "score", 
              "bias", 
              "#", 
              "of", 
              "c-Evalue", 
              "i-Evalue", 
              "score", 
              "bias",
              "from",
              "to",
              "from",
              "to",
              "from",
              "to",
              "acc",
              "description of target")


#remove hashed lines 
HMM_results <- sapply(HMM_results, simplify = FALSE, USE.NAMES = TRUE, FUN = function(i) {
  df <- i[!grepl("#", i$V1),]
  colnames(df)<- col_names
  return(df)
})


#parse by e-value
#get only hits with c-Evalue < e-3
HMM_results_high_confidence<- sapply(HMM_results, simplify = FALSE, USE.NAMES = TRUE, FUN = function(i) {
  df <- i[as.numeric(i$'E-value') < 1E-10,]
  colnames(df)<- col_names
  return(df)
})



#split the dataframes 
split_the_data_frames<- function(df){
  #split dfs
  split_dfs<- split(df, f = df$`query name`)
}


#run function to split all data frames in the list
split_dfs2<- lapply(HMM_results_high_confidence, function(x) split_the_data_frames(x))


#function to identify incongruencies:
#targest with hits in more than one class, are reassigned to "Unknown"
rename_incongruencies<- function(lsted){
  out_list <- lsted
  for (i in 1:length(lsted)){
    for (j in 1:length(lsted[[i]])){
      if (nrow(unique(split_dfs2[[i]][[j]][,1])) >1) {
        lsted[[i]][[j]][,1] <- "Unknown"
        out_list[[i]][[j]][,1] <- lsted[[i]][[j]][,1]
      }
    }
  }
  return(out_list)
}

#instanciate list
renamed_list<- list()
#run the function
renamed_list<- rename_incongruencies(lsted= split_dfs2) 


#function to merge dataframes back into one df per species 
merge_dfs<- function(lsted){
  ans<- list()
  out_list <- list()
  for (i in 1:length(lsted)){
    #for (j in 1:length(lsted[[i]])){
    out_list<- list(rbindlist(lsted[[i]]))
    ans <- append(ans, out_list)
  }
  return(ans)
}

#run function 
merged_dfs<- merge_dfs(lsted = renamed_list)

#shrink to only unique lines for each query
one_unique_represent <-lapply(merged_dfs, function(i) i[!duplicated(i$'query name'), ])

##seperate Suillus from Other
#cobine all data frames
one_big_df<- data.frame(bind_rows(one_unique_represent, .id = "Species"))
#split based on Sui vs. Other
Suillus<- one_big_df[grep("*Sui*", one_big_df$query.name), ]
Other<- one_big_df[grep("*Sui*", one_big_df$query.name, invert = TRUE), ]

#call table based on the Species col and the target name col
Suillus_table<- table(Suillus$Species, Suillus$target.name)
Other_table<- table(Other$Species, Other$target.name)

Suillus_table2<- as.data.frame(Suillus_table)
Other_table2<- as.data.frame(Other_table)
###make figures

#set colnames
colnames(Suillus_table2)<- c("Species", "Class", "Count")
colnames(Other_table2)<- c("Species", "Class", "Count")

#add a group collumn 
Suillus_table2$Group<- "S"
Other_table2$Group<- "O"

##for S vs O:
#bind the groups back together 
input_df1<- rbind(Suillus_table2, Other_table2)
#remove "GPCR_class" for sorting
input_df2 <- data.frame(lapply(input_df1, function(x) {
  gsub("GPCR_class", "", x)
}))
#order by class 
input_df3<- input_df2[order(as.numeric(as.character(input_df2$Class))), ]
#change back to chr. 
input_df3$Count<- as.integer(as.character(input_df3$Count))
#set factors

#input_df3$Class2 <- factor(input_df3$Class, as.character(input_df3$Class))
input_df3$Class2 <- factor(input_df3$Class, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, "Unknown"))

##make grouped strip chart
p<- ggplot(input_df3, aes(x=Class2, y=Count, fill=Group, color = Group)) + 
  geom_point(position=position_jitterdodge(jitter.width = .15), aes(shape=Group, color=Group), size=2, alpha=0.5)+
  scale_shape_manual(values=c(19, 19, 19))+
  scale_color_manual(values=c("#5676A1", "#B09136", "#B36757"))+
  theme(legend.position="top") +
  stat_summary(fun.y = mean, geom = "point", aes(group = interaction(Class2, Group)),
               shape = 95, size = 10, show.legend = F, alpha=1)+
  stat_summary(geom = "point",fun.y = "mean",col = "black", size = 5, shape = 2,fill = "black", alpha=.8, show.legend = F,) +
  theme_minimal(base_size = 7) +
  scale_x_discrete(labels=c("Class 1","Class 2","Class 3","Class 4","Class 5", "Class 8", "Class 9", "Class 10","Class 11","Class 12","Class 13","Class 14","Unknown"))+
  scale_y_discrete(limits=c(1,3,5,7,9,11,13,15))
dev.new()
p + coord_equal(ratio = .3) +
  theme(axis.title.x = element_blank()) +
  labs(y = "n GPCRs")




###you are here
#overall - ave n GPCRs per species. (Suillus vs. Other) 
#function
split_the_data_frames_by_species<- function(df){split(df, f = df$Species)
}
#run function
split_dfs_Suillus<- split_the_data_frames_by_species(df = Suillus_table2)
split_dfs_Other<- split_the_data_frames_by_species(df = Other_table2)

#get averages for each species 
suillus_GPCRs<- lapply(split_dfs_Suillus, function(x) sum(x$Count))
suillus_GPCRs2<- round(sum(as.numeric(suillus_GPCRs)) / length(suillus_GPCRs), 0)
other_GPCRs<- lapply(split_dfs_Other, function(x) sum(x$Count))
other_GPCRs2<- round(sum(as.numeric(other_GPCRs)) / length(other_GPCRs), 0)


all_GPCRs<- data.frame(rbind(Suillus= suillus_GPCRs2, Other= other_GPCRs2))
colnames(all_GPCRs)<- "Count"
all_GPCRs$Group<- row.names(all_GPCRs)

#reorder
all_GPCRs$Group <- factor(all_GPCRs$Group, c("Suillus", "Other"))

#add sd
suillus_sd<- lapply(split_dfs_Suillus, function(x) sd(x$Count))
suillus_sd2<- sum(as.numeric(suillus_sd) / length(suillus_sd))
other_sd<- lapply(split_dfs_Other, function(x) sd(x$Count))
other_sd2<- sum(as.numeric(other_sd) / length(other_sd))

all_GPCRs$sd<- rbind(suillus_sd2, other_sd2)


#plot
p<-ggplot(all_GPCRs, aes(x=Group, y=Count, fill=Group)) +
  geom_bar(stat="identity") +
  theme_minimal()
p + scale_fill_manual(values=c("#5676A1", "#B09136", "#B36757")) +
  theme(axis.title.x = element_blank()) +
  labs(y = "ave. n GPCRs") +
  geom_errorbar(aes(ymin=Count-sd, ymax=Count+sd), width=.2,
                position=position_dodge(.9)) 












##for within S compar.
#create groups for each host association 
White<- c("Suiame1", 
          "Suipla1", 
          "Suipic1") 

Red<- c("Suibov1", 
        "Suibr2", 
        "Suicot1", 
        "Suidec1", 
        "Suifus1", 
        "Suihi1", 
        "Suilu4", 
        "Suiocc1", 
        "Suisu1", 
        "Suitom1", 
        "Suivar1")

Larch<- c("Suiamp1", 
          "Suicli1", 
          "Suipal1")


#seperate data frames based on host association 
White_counts<- Suillus[grepl(paste(White, collapse="|"), Suillus$query.name), ]
Red_counts<- Suillus[grepl(paste(Red, collapse="|"), Suillus$query.name), ]
Larch_counts<- Suillus[grepl(paste(Larch, collapse="|"), Suillus$query.name), ]

##plot
#call table
Suillus_table_white<- data.frame(table(White_counts$Species, White_counts$target.name))
Suillus_table_red<- data.frame(table(Red_counts$Species, Red_counts$target.name))
Suillus_table_larch<- data.frame(table(Larch_counts$Species, Larch_counts$target.name))
colnames(Suillus_table_white)<- c("Species", "Class", "Count")
colnames(Suillus_table_red)<- c("Species", "Class", "Count")
colnames(Suillus_table_larch)<- c("Species", "Class", "Count")

#add host group designator 
Suillus_table_white$host_group<- "W"
Suillus_table_red$host_group<- "R"
Suillus_table_larch$host_group<- "L"

#bind them back together
Suillus_df<- rbind(Suillus_table_white, Suillus_table_red, Suillus_table_larch)
#remove any zero values 
Suillus_df_no_zeros<- Suillus_df[Suillus_df$Count != 0,]

#remove "GPCR_class" for sorting
input_df2_S <- data.frame(lapply(Suillus_df_no_zeros, function(x) {
  gsub("GPCR_class", "", x)
}))
#order by class 
input_df3_S<- input_df2_S[order(as.numeric(as.character(input_df2_S$Class))), ]
#change back to chr. 
input_df3_S$Count<- as.integer(as.character(input_df3_S$Count))
#set factors

input_df3_S$Class2 <- factor(input_df3_S$Class, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, "Unknown"))



#colors: pink: B36757, yellow: B09136, blue: 5676A1
#get means 

##make grouped strip chart
p<- ggplot(input_df3_S, aes(x=Class2, y=Count, fill=host_group, color = host_group)) + 
  geom_point(position=position_jitterdodge(jitter.width = .15), aes(shape=host_group, color=host_group), size=2, alpha=0.5)+
  scale_shape_manual(values=c(19, 19, 19))+
  scale_color_manual(values=c("#5676A1", "#B09136", "#B36757"))+
  theme(legend.position="top") +
  stat_summary(fun.y = mean, geom = "point", aes(group = interaction(Class2, host_group)),
               shape = 95, size = 10, show.legend = F, alpha=1)+
  stat_summary(geom = "point",fun.y = "mean",col = "black", size = 5, shape = 2,fill = "black", alpha=.8, show.legend = F,) +
  theme_minimal(base_size = 8) +
  scale_x_discrete(labels=c("Class 1","Class 2","Class 3","Class 4","Class 8","Class 10","Class 11","Class 12","Class 13","Class 14","Unknown"))+
  scale_y_discrete(limits=c(1,3,5,7,9,11))
#scale_x_discrete(expand=c(.1, 0))+
#scale_y_discrete(expand=c(.001, .1))
dev.new()
p + coord_equal(ratio = .3) +
  theme(axis.title.x = element_blank()) +
  labs(y = "n GPCRs")



#overall - ave n GPCRs per species. 
#function
split_the_data_frames_by_species<- function(df){split(df, f = df$Species)
}
#run function
split_dfs_white<- split_the_data_frames_by_species(df = Suillus_table_white)
split_dfs_red<- split_the_data_frames_by_species(df = Suillus_table_red)
split_dfs_larch<- split_the_data_frames_by_species(df = Suillus_table_larch)

#get averages for each species 
white_GPCRs<- lapply(split_dfs_white, function(x) sum(x$Count))
white_GPCRs2<- round(sum(as.numeric(white_GPCRs)) / length(white_GPCRs), 0)
red_GPCRs<- lapply(split_dfs_red, function(x) sum(x$Count))
red_GPCRs2<- round(sum(as.numeric(red_GPCRs)) / length(red_GPCRs), 0)
larch_GPCRs<- lapply(split_dfs_larch, function(x) sum(x$Count))
larch_GPCRs2<- round(sum(as.numeric(larch_GPCRs)) / length(larch_GPCRs), 0)

S_GPCRs<- data.frame(rbind( Red= red_GPCRs2, White= white_GPCRs2, Larch=larch_GPCRs2))
colnames(S_GPCRs)<- "Count"
S_GPCRs$Group<- row.names(S_GPCRs)

#reorder
S_GPCRs$Group <- factor(S_GPCRs$Group, c("Red", "White", "Larch"))

#add sd
white_sd<- lapply(split_dfs_white, function(x) sd(x$Count))
white_sd2<- sum(as.numeric(white_sd) / length(white_sd))
red_sd<- lapply(split_dfs_red, function(x) sd(x$Count))
red_sd2<- sum(as.numeric(red_sd) / length(red_sd))
larch_sd<- lapply(split_dfs_larch, function(x) sd(x$Count))
larch_sd2<- sum(as.numeric(larch_sd) / length(larch_sd))
S_GPCRs$sd<- rbind(red_sd2, white_sd2, larch_sd2)


#plot
p<-ggplot(S_GPCRs, aes(x=Group, y=Count, fill=Group)) +
  geom_bar(stat="identity") +
  theme_minimal()
p + scale_fill_manual(values=c("#5676A1", "#B09136", "#B36757")) +
  theme(axis.title.x = element_blank()) +
  labs(y = "ave. n GPCRs") +
  geom_errorbar(aes(ymin=Count-sd, ymax=Count+sd), width=.2,
                position=position_dodge(.9)) 
