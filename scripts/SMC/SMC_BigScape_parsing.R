
#This script looks at the distribution of BigScape identified Terpene and NRPS secondary metabolite clusters
setwd("~/Desktop/BIGSCAPE_results_feb_2020_pub_version/")

#load packages
library("seqinr")
library("tidyr")
library("UpSetR")

#read in input files: (these were created using Jason's Parsing script that creates counts from the 'Orthogroups' output file from orthofinder)
#Terpenes_Suillus<- read.table("absence_presence_terp_suillus_only.tab", header = TRUE)

#new version
Terpenes_Suillus<- read.table("absence_presence_Terpenes.tab", header = TRUE)

#white list
White_list<- c("S_americanus", 
               "S_pictus", 
               "S_placedus", 
               "S_plorans",
               "S_subluteus",
               "S_discolor")

Red_list<- c("S_bovinus", 
             "S_brevipes", 
             "S_cothurnatus", 
             "S_decipiens", 
             "S_fuscotomentosus", 
             "S_hirtellus", 
             "S_luteus", 
             "S_occidentalis", 
             "S_subalutaceus", 
             "S_tomentosus", 
             "S_variegatus")

Larch_list<- c("S_ampliporus", 
               "S_paluster", 
               "S_clintonianus")



#if it's not zero in column in the<>list, retun that collumn  
White_df<- Terpenes_Suillus[Terpenes_Suillus$ACC %in% White_list, ]
row.names(White_df)<- White_df$ACC
White_df<- White_df[,2:length(White_df)]

Red_df<- Terpenes_Suillus[Terpenes_Suillus$ACC %in% Red_list, ]
row.names(Red_df)<- Red_df$ACC
Red_df<- Red_df[,2:length(Red_df)]

Larch_df<- Terpenes_Suillus[Terpenes_Suillus$ACC %in% Larch_list, ]
row.names(Larch_df)<- Larch_df$ACC
Larch_df<- Larch_df[,2:length(Larch_df)]

#change original
row.names(Terpenes_Suillus)<- Terpenes_Suillus$ACC
Terpenes_Suillus<- Terpenes_Suillus[,2:length(Terpenes_Suillus)]

#return the collumn in the whole data set if the row is not zero in the <>df
based_on_white<- Terpenes_Suillus[,colSums(White_df != 0) > 0]
based_on_red<- Terpenes_Suillus[,colSums(Red_df != 0) > 0]
based_on_larch<- Terpenes_Suillus[,colSums(Larch_df != 0) > 0]

#remove host group from the subsest dfs
lacking_white<- based_on_white[! rownames(based_on_white) %in% White_list, ]
lacking_red<- based_on_red[! rownames(based_on_red) %in% Red_list, ]
lacking_larch<- based_on_larch[! rownames(based_on_larch) %in% Larch_list, ]

#get host only
#white_only<- based_on_white[rownames(based_on_white) %in% White_list, ]
#red_only<- based_on_red[rownames(based_on_red) %in% Red_list, ]
#larch_only<- based_on_larch[rownames(based_on_larch) %in% Larch_list, ]

#return all collumns that are == 0 when the host group is taken out
unique_to_white<- lacking_white[, apply(lacking_white == 0, 2, all)]
unique_to_red<- lacking_red[, apply(lacking_red == 0, 2, all)]
unique_to_larch<- lacking_larch[, apply(lacking_larch == 0, 2, all)]

dim(unique_to_white)
dim(unique_to_red)
dim(unique_to_larch)
dim(Terpenes_Suillus)

#are there any Terpene clusters that are in ALL the representatives of that host but no others?

#if it's not zero for ALL the columns in the<>list, retun that collumn  
#look to make sure they all have one
#White_df[,colSums(White_df) >= nrow(White_df)]
#Red_df[,colSums(Red_df) >= nrow(Red_df)]
#Larch_df[,colSums(Larch_df) >= nrow(Larch_df)]

#looks good, subset
new_white<- White_df[,colSums(White_df) >= nrow(White_df)]
new_red<-Red_df[,colSums(Red_df) >= nrow(Red_df)]
new_larch<-Larch_df[,colSums(Larch_df) >= nrow(Larch_df)]

#return all rows in the whole data set that are in the col names in new<>
based_on_white2<- Terpenes_Suillus[,colnames(new_white)]
based_on_red2<- Terpenes_Suillus[,colnames(new_red)]
based_on_larch2<- Terpenes_Suillus[,colnames(new_larch)]

#are any clusters in all of one host group and none of the others?
#(skipping- we can see from the above output that there are not)



###Make upset plot
#change to numeric
Terpenes_Suillus_numeric<- as.data.frame(lapply(Terpenes_Suillus, as.numeric))
class(Terpenes_Suillus_numeric[1,1])

#change all 2's and 3's to 1 to make binary presence / absence
Terpenes_Suillus_binary<- as.data.frame(lapply(Terpenes_Suillus_numeric, gsub, pattern = 2, replacement = 1, fixed = TRUE))
Terpenes_Suillus_binary_2<- as.data.frame(lapply(Terpenes_Suillus_binary, gsub, pattern = 3, replacement = 1, fixed = TRUE))
row.names(Terpenes_Suillus_binary_2)<- rownames(Terpenes_Suillus)

sum(Terpenes_Suillus_numeric == 2)
sum(Terpenes_Suillus_binary == 2)
sum(Terpenes_Suillus_binary == 3)
sum(Terpenes_Suillus_binary_2 == 3)

#split by host association
White_df<- Terpenes_Suillus_binary_2[rownames(Terpenes_Suillus_binary_2) %in% White_list, ]
Red_df<- Terpenes_Suillus_binary_2[rownames(Terpenes_Suillus_binary_2) %in% Red_list, ]
Larch_df<- Terpenes_Suillus_binary_2[rownames(Terpenes_Suillus_binary_2) %in% Larch_list, ]

#chagne back to numeric
White_df<- as.data.frame(lapply(White_df, as.numeric))
Red_df<- as.data.frame(lapply(Red_df, as.numeric))
Larch_df<- as.data.frame(lapply(Larch_df, as.numeric))

#get averages
White<- data.frame(White= round(colMeans(White_df), 0))
Red<- data.frame(Red= round(colMeans(Red_df), 0))
Larch<- data.frame(Larch= round(colMeans(Larch_df), 0))

input_df<- cbind(White, Red, Larch)

#make upset plot of terpenes
upset(input_df, sets = c("Larch", "White", "Red"), mb.ratio = c(0.55, 0.45), order.by = "degree",
      mainbar.y.label = "Intersection Size - Terpenes",
      keep.order = TRUE)
# no particular trends 

###do the same plot for NRPS
#NRPS_Suillus_init<- read.table("absence_presence_NRPS_Suillus_only.tab", header = TRUE)
NRPS_Suillus_init<- read.table("absence_presence_NRPS.tab", header = TRUE)
NRPS_Suillus<- NRPS_Suillus_init[, 2:ncol(NRPS_Suillus_init)]
row.names(NRPS_Suillus)<- NRPS_Suillus_init$ACC

#change to numeric
NRPS_Suillus_numeric<- as.data.frame(lapply(NRPS_Suillus, as.numeric))
class(NRPS_Suillus_numeric[1,1])

#change all 2's and 3's to 1 to make binary presence / absence
NRPS_Suillus_binary<- as.data.frame(lapply(NRPS_Suillus_numeric, gsub, pattern = 2, replacement = 1, fixed = TRUE))
NRPS_Suillus_binary_2<- as.data.frame(lapply(NRPS_Suillus_binary, gsub, pattern = 3, replacement = 1, fixed = TRUE))
NRPS_Suillus_binary_3<- as.data.frame(lapply(NRPS_Suillus_binary_2, gsub, pattern = 4, replacement = 1, fixed = TRUE))
NRPS_Suillus_binary_4<- as.data.frame(lapply(NRPS_Suillus_binary_3, gsub, pattern = 5, replacement = 1, fixed = TRUE))
NRPS_Suillus_binary_5<- as.data.frame(lapply(NRPS_Suillus_binary_4, gsub, pattern = 6, replacement = 1, fixed = TRUE))


row.names(NRPS_Suillus_binary_5)<- rownames(NRPS_Suillus)

sum(NRPS_Suillus_numeric == 2)
sum(NRPS_Suillus_binary == 2)
sum(NRPS_Suillus_binary == 3)
sum(NRPS_Suillus_binary_2 == 3)
sum(NRPS_Suillus_binary_3 == 4)
sum(NRPS_Suillus_numeric == 5)
sum(NRPS_Suillus_binary_4 == 5)
sum(NRPS_Suillus_numeric == 6)
sum(NRPS_Suillus_binary_5 == 6)

#split by host association
White_df<- NRPS_Suillus_binary_5[rownames(NRPS_Suillus_binary_5) %in% White_list, ]
Red_df<- NRPS_Suillus_binary_5[rownames(NRPS_Suillus_binary_5) %in% Red_list, ]
Larch_df<- NRPS_Suillus_binary_5[rownames(NRPS_Suillus_binary_5) %in% Larch_list, ]

#chagne back to numeric
White_df_2<- as.data.frame(lapply(White_df, as.numeric))
Red_df_2<- as.data.frame(lapply(Red_df, as.numeric))
Larch_df_2<- as.data.frame(lapply(Larch_df, as.numeric))

#get averages 
White<- data.frame(White= round(colMeans(White_df_2), 0))
Red<- data.frame(Red= round(colMeans(Red_df_2), 0))
Larch<- data.frame(Larch= round(colMeans(Larch_df_2), 0))

input_df2<- cbind(White, Red, Larch)

#make upset plot of NRPS
upset(input_df2, sets = c("Larch", "White", "Red"), mb.ratio = c(0.55, 0.45), order.by = "degree",
      mainbar.y.label = "Intersection Size - NRPS",
      keep.order = TRUE)

