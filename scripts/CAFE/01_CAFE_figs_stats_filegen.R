#processing script for CAFE results files

#setwd
setwd("~/Desktop/CAFE_output_Feb_2020")

#load libraries
library("coin")
library("multcompView")
library("multcomp")
library("car")

#to load coin
#if(!require(coin)){install.packages("coin")}
#if(!require(FSA)){install.packages("FSA")}
#if(!require(rcompanion)){install.packages("rcompanion")}
#if(!require(multcompView)){install.packages("multcompView")}
options(stringsAsFactors = FALSE)

#read in files
summary<- read.csv("summary_NEW_r0_pub.csv", header = TRUE)


#split groups 
White_list<- c("Suiame1", 
               "Suipic1", 
               "Suipla1", 
               "Suiplo1",
               "Suisubl1",
               "Suidis1")

Red_list<- c("Suibov1", 
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

Larch_list<- c("Suiamp1", 
               "Suipal1", 
               "Suicli1")

White_df<- summary[summary$Species %in% White_list,]
Red_df<- summary[summary$Species %in% Red_list,]
Larch_df<- summary[summary$Species %in% Larch_list,]

Suillus_df<- summary[grep("Sui", summary$Species),]
Other_df<- summary[grep("Sui", summary$Species, invert = TRUE),]
suspects<- suspects<- c("Gyrli1", "Rhivul1", "Rhivi1", "Rhives1", "Rhisa1")
Other_df_wo<- Other_df[!Other_df$Species %in% suspects,]


#add designators and recombine 
White_df$group<- "W"
Red_df$group<- "R"
Larch_df$group<- "L"

Suillus_df$group<- "S"
Other_df$group<- "O"
Other_df_wo$group<- "O"

Suillus_df2<- rbind(White_df, Red_df, Larch_df)
Other_df2<- rbind(Suillus_df, Other_df)
Other_df2_wo<- rbind(Suillus_df, Other_df_wo)


##S vs. O
#t_test for expansions


#test for normality
shapiro.test(as.numeric(Other_df2$sig_Expanded_fams) [Other_df2$group == "S"])
#Suillus not normal
shapiro.test(as.numeric(log(Other_df2$sig_Expanded_fams)) [Other_df2$group == "S"])
#still not normal with log transform
powerTransform(Other_df2$sig_Expanded_fams)
#recomended -0.9496181
shapiro.test(as.numeric((Other_df2$sig_Expanded_fams)**-1) [Other_df2$group == "S"])

shapiro.test(as.numeric(Other_df2$sig_Expanded_fams) [Other_df2$group == "O"])
#Other is normal, but transform to match 
shapiro.test(as.numeric((Other_df2$sig_Expanded_fams)**-1) [Other_df2$group == "O"])

Other_df2$sig_Expanded_fams_transform<- Other_df2$sig_Expanded_fams**-1


#test for equal variance 
var.test(as.numeric(Other_df2$sig_Expanded_fams_transform [Other_df2$group == "S"]), 
         (as.numeric(Other_df2$sig_Expanded_fams_transform) [Other_df2$group == "O"]))
#fail to reject, variance is equalish


t.test(as.numeric(sig_Expanded_fams_transform) ~ as.character(group), data = Other_df2, var.equal = TRUE)
O_E<- mean(Other_df$sig_Expanded_fams)
S_E<- mean(Suillus_df$sig_Expanded_fams)
#suillus has significantly more gene family expansions

##excluding High and Moderare 
#test for normality
shapiro.test(as.numeric(Other_df2_wo$sig_Expanded_fams) [Other_df2_wo$group == "S"])
#Suillus not normal
shapiro.test(as.numeric(log(Other_df2_wo$sig_Expanded_fams)) [Other_df2_wo$group == "S"])
#still not normal with log transform
powerTransform(Other_df2_wo$sig_Expanded_fams)
#recomended -0.9496181
shapiro.test(as.numeric((Other_df2_wo$sig_Expanded_fams)**-1) [Other_df2_wo$group == "S"])

shapiro.test(as.numeric(Other_df2_wo$sig_Expanded_fams) [Other_df2_wo$group == "O"])
#Other is normal, but transform to match 
shapiro.test(as.numeric((Other_df2_wo$sig_Expanded_fams)**-1) [Other_df2_wo$group == "O"])

Other_df2_wo$sig_Expanded_fams_transform<- Other_df2_wo$sig_Expanded_fams**-1


#test for equal variance 
var.test(as.numeric(Other_df2_wo$sig_Expanded_fams_transform [Other_df2_wo$group == "S"]), 
         (as.numeric(Other_df2_wo$sig_Expanded_fams_transform) [Other_df2_wo$group == "O"]))
#fail to reject, variance is equalish


t.test(as.numeric(sig_Expanded_fams_transform) ~ as.character(group), data = Other_df2_wo, var.equal = TRUE)
O_E<- mean(Other_df_wo$sig_Expanded_fams)
S_E<- mean(Suillus_df$sig_Expanded_fams)
#suillus has significantly more gene family expansions


#t_test for CONTRACTIONS
#test for normality
shapiro.test(as.numeric(Other_df2$sig_Contracted_fams) [Other_df2$group == "S"])
#not normal
shapiro.test(as.numeric(Other_df2$sig_Contracted_fams) [Other_df2$group == "O"])
#not normal - need to transform

#test for equal variance 
var.test(as.numeric(Other_df2$sig_Contracted_fams [Other_df2$group == "S"]), 
         (as.numeric(Other_df2$sig_Contracted_fams) [Other_df2$group == "O"]))
#variance not equal. 

#transform 
Other_df2$log_sig_Contracted_fams<- log(Other_df2$sig_Contracted_fams+1)


#did that help?
#test for normality
shapiro.test(as.numeric(Other_df2$log_sig_Contracted_fams) [Other_df2$group == "S"])
shapiro.test(as.numeric(Other_df2$log_sig_Contracted_fams) [Other_df2$group == "O"])
#yep, normal now. 
#retest var.
var.test(as.numeric(Other_df2$log_sig_Contracted_fams [Other_df2$group == "S"]), 
         (as.numeric(Other_df2$log_sig_Contracted_fams) [Other_df2$group == "O"]))
#variance not equal. 

t.test(as.numeric(log_sig_Contracted_fams) ~ as.character(group), data = Other_df2, var = FALSE)
O_C<- mean(Other_df$sig_Contracted_fams)
S_C<- mean(Suillus_df$sig_Contracted_fams)
#suillus has significantly more gene family contractions too

#removing high and moderate from the dataset
#t_test for CONTRACTIONS
#test for normality
shapiro.test(as.numeric(Other_df2_wo$sig_Contracted_fams) [Other_df2_wo$group == "S"])
#not normal
shapiro.test(as.numeric(Other_df2_wo$sig_Contracted_fams) [Other_df2_wo$group == "O"])
#not normal - need to transform

#test for equal variance 
var.test(as.numeric(Other_df2_wo$sig_Contracted_fams [Other_df2_wo$group == "S"]), 
         (as.numeric(Other_df2_wo$sig_Contracted_fams) [Other_df2_wo$group == "O"]))
#variance not equal. 

#transform 
Other_df2_wo$log_sig_Contracted_fams<- log(Other_df2_wo$sig_Contracted_fams+1)


#did that help?
#test for normality
shapiro.test(as.numeric(Other_df2_wo$log_sig_Contracted_fams) [Other_df2_wo$group == "S"])
shapiro.test(as.numeric(Other_df2_wo$log_sig_Contracted_fams) [Other_df2_wo$group == "O"])
#yep, normal-ish now. 
#retest var.
var.test(as.numeric(Other_df2_wo$log_sig_Contracted_fams [Other_df2_wo$group == "S"]), 
         (as.numeric(Other_df2_wo$log_sig_Contracted_fams) [Other_df2_wo$group == "O"]))
#variance not equal. 

t.test(as.numeric(log_sig_Contracted_fams) ~ as.character(group), data = Other_df2_wo, var = FALSE)
O_C<- mean(Other_df_wo$sig_Contracted_fams)
S_C<- mean(Suillus_df$sig_Contracted_fams)

#make figure:
## S vs. O
#make new df with negatives 
sig_Contracted_fams_SO<- cbind(as.numeric(-Other_df2$sig_Contracted_fams), Other_df2$group)
sig_Expanded_fams_SO<- cbind(as.numeric(Other_df2$sig_Expanded_fams), Other_df2$group)
df_for_graph_S_vs_O<- rbind(sig_Contracted_fams_SO, sig_Expanded_fams_SO)
colnames(df_for_graph_S_vs_O)<- c("val", "group")
df_for_graph_S_vs_O<- as.data.frame(df_for_graph_S_vs_O)
df_for_graph_S_vs_O$group = factor(df_for_graph_S_vs_O$group,c("S","O"))
x<- stripchart(as.numeric(val) ~ group, data = df_for_graph_S_vs_O,
               ylab="n gene families contracted",
               method="jitter",
               col=c("#40595266","#9A8F5A66"),
               pch=19, 
               ylim= c(-120, 120),
               vertical=TRUE
)
abline(a=1, b= 0, col = "black", lty = 3)

#add means
segments(.85,S_E, 1.15,S_E, col = "#405952", lwd = 2)
segments(1.85,O_E, 2.15, O_E, col = "#9A8F5A", lwd = 2)

segments(.85,-S_C, 1.15,-S_C, col = "#405952", lwd = 2)
segments(1.85,-O_C, 2.15, -O_C, col = "#9A8F5A", lwd = 2)





#Suillus has greater expansions and contractions


### W vs. R vs. L

#vis normality and variance
#build model
#expanded fams
sig_Expanded_fams_model <- lm(Suillus_df2$sig_Expanded_fams ~ Suillus_df2$group)
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(sig_Expanded_fams_model) #data is non normal, with non-constant variance 
par(mfrow = c(1,1)) #back to one image per page
boxCox(sig_Expanded_fams_model) #what does boxcox suggest?  normalish, bit non constant var.
powerTransform(sig_Expanded_fams_model) 
Suillus_df2$log_sig_Expanded_fams<- log(Suillus_df2$sig_Expanded_fams)
Suillus_df2$log_sig_Expanded_neghalf<- Suillus_df2$sig_Expanded_fams^(-.5)
#did it help? 
log_sig_Expanded_fams_model <- lm(Suillus_df2$log_sig_Expanded_fams ~ Suillus_df2$group)
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(log_sig_Expanded_fams_model) #data is non normal, with non-constant variance 
par(mfrow = c(1,1)) #back to one image per page

#can not meet variance assumptions regardless of transformation - randomization test.
#randomization test make df
r_GS<- data.frame(Red_df$sig_Expanded_fams, rep("Red"))
colnames(r_GS)<- c("GS", "group")
w_GS<- data.frame(White_df$sig_Expanded_fams, rep("White"))
colnames(w_GS)<- c("GS", "group")
l_GS<- data.frame(Larch_df$sig_Expanded_fams, rep("Larch"))
colnames(l_GS)<- c("GS", "group")
#bind them
GS_df<- rbind(r_GS, w_GS, l_GS)
#set factors
#set order of boxes
GS_df$group <- factor(GS_df$group,levels = c('Red','White', 'Larch'),ordered = TRUE)

#run independence test (this is like an ANOVA)
independence_test(GS ~ group, 
                  data = GS_df)

#not significant! 
#Z = -1.3355, p-value = 0.1817


##contracted fams
#randomization test make df
r_GS<- data.frame(Red_df$Contracted_fams, rep("Red"))
colnames(r_GS)<- c("GS", "group")
w_GS<- data.frame(White_df$Contracted_fams, rep("White"))
colnames(w_GS)<- c("GS", "group")
l_GS<- data.frame(Larch_df$Contracted_fams, rep("Larch"))
colnames(l_GS)<- c("GS", "group")
#bind them
GS_df<- rbind(r_GS, w_GS, l_GS)
#set factors
#set order of boxes
GS_df$group <- factor(GS_df$group,levels = c('Red','White', 'Larch'),ordered = TRUE)

#run independence test (this is like an ANOVA)
independence_test(GS ~ group, 
                  data = GS_df, distribution = approximate(nresample = 10000))



#not significant
mean(White_df$sig_Expanded_fams)
mean(Red_df$sig_Expanded_fams)
mean(Larch_df$sig_Expanded_fams)

summary(aov(sig_Contracted_fams ~ group, data = Suillus_df2))
#significant - but not relevent, does not meet assumptions

mean(White_df$sig_Contracted_fams)
mean(Red_df$sig_Contracted_fams)
mean(Larch_df$sig_Contracted_fams)

#make figure for R vs W vs L
#make new df with negatives 
sig_Contracted_fams_RWL<- cbind(as.numeric(-Suillus_df2$sig_Contracted_fams), Suillus_df2$group)
sig_Expanded_fams_RWL<- cbind(as.numeric(Suillus_df2$sig_Expanded_fams), Suillus_df2$group)
df_for_graph_RWL<- rbind(sig_Contracted_fams_RWL, sig_Expanded_fams_RWL)
colnames(df_for_graph_RWL)<- c("val", "group")
df_for_graph_RWL<- as.data.frame(df_for_graph_RWL)
df_for_graph_RWL$group = factor(df_for_graph_RWL$group,c("R","W", "L"))

x<- stripchart(as.numeric(val) ~ group, data = df_for_graph_RWL,
               ylab="n gene families contracted",
               method="jitter",
               col=c("#5676A166","#B0913666", "#B3675766"),
               pch=19, 
               ylim= c(-120, 120),
               vertical=TRUE
)
abline(a=1, b= 0, col = "black", lty = 3)
#add mean lines
segments(.75,mean(Red_df$sig_Expanded_fams), 1.25,mean(Red_df$sig_Expanded_fams), col = "#5676A1", lwd = 2)
segments(1.75,mean(White_df$sig_Expanded_fams), 2.25,mean(White_df$sig_Expanded_fams), col = "#B09136", lwd = 2)
segments(2.75,mean(Larch_df$sig_Expanded_fams), 3.25,mean(Larch_df$sig_Expanded_fams), col = "#B36757", lwd = 2)

segments(.75,mean(-Red_df$sig_Contracted_fams), 1.25,mean(-Red_df$sig_Contracted_fams), col = "#5676A1", lwd = 2)
segments(1.75,mean(-White_df$sig_Contracted_fams), 2.25,mean(-White_df$sig_Contracted_fams), col = "#B09136", lwd = 2)
segments(2.75,mean(-Larch_df$sig_Contracted_fams), 3.25,mean(-Larch_df$sig_Contracted_fams), col = "#B36757", lwd = 2)



###specific fams analysis 
summary_fams_All<- read.table(text = gsub(",", "\t", readLines("summary_NEW_r0_fams.txt")), header = FALSE, sep = "\t", check.names = FALSE, fill = TRUE)
summary_fams_All<- summary_fams_All[-1,]

#remove inner node information: 
list<- paste0("<", rep(1:100), ">:")

#all_leaves_Suillus<- summary_fams_Suillus[!(summary_fams_Suillus$V1) %in% list,]
all_leaves_All<- summary_fams_All[!(summary_fams_All$V1) %in% list,]

#group from all - this is the group you're actually interested in given the first set of significance tests
Suillus_only<- data.frame(all_leaves_All[grep("Sui", all_leaves_All$V1),])
Suillus_none<- data.frame(all_leaves_All[grep("Sui", all_leaves_All$V1, invert = TRUE),])

White<- data.frame(Suillus_only[ grepl(paste(White_list,  collapse="|"), Suillus_only$V1),])
Red<- data.frame(Suillus_only[ grepl(paste(Red_list,  collapse="|"), Suillus_only$V1),])
Larch<- data.frame(Suillus_only[ grepl(paste(Larch_list,  collapse="|"), Suillus_only$V1),])


#seperate expansions and contractions
Suillus_only_C <- list()
for (i in 1:nrow(Suillus_only)){
  Suillus_only_C[[i]]<- grep("-",Suillus_only[i,], value = TRUE, fixed = TRUE)
}

Suillus_only_E <- list()
for (i in 1:nrow(Suillus_only)){
  Suillus_only_E[[i]]<- grep("+",Suillus_only[i,], value = TRUE, fixed = TRUE)
}


Suillus_none_C <- list()
for (i in 1:nrow(Suillus_none)){
  Suillus_none_C[[i]]<- grep("-",Suillus_none[i,], value = TRUE, fixed = TRUE)
}

Suillus_none_E <- list()
for (i in 1:nrow(Suillus_none)){
  Suillus_none_E[[i]]<- grep("+",Suillus_none[i,], value = TRUE, fixed = TRUE)
}

Red_E <- list()
for (i in 1:nrow(Red)){
  Red_E[[i]]<- grep("+",Red[i,], value = TRUE, fixed = TRUE)
}
White_E <- list()
for (i in 1:nrow(White)){
  White_E[[i]]<- grep("+",White[i,], value = TRUE, fixed = TRUE)
}
Larch_E <- list()
for (i in 1:nrow(Larch)){
  Larch_E[[i]]<- grep("+",Larch[i,], value = TRUE, fixed = TRUE)
}

Red_C <- list()
for (i in 1:nrow(Red)){
  Red_C[[i]]<- grep("-",Red[i,], value = TRUE, fixed = TRUE)
}
White_C <- list()
for (i in 1:nrow(White)){
  White_C[[i]]<- grep("-",White[i,], value = TRUE, fixed = TRUE)
}
Larch_C <- list()
for (i in 1:nrow(Larch)){
  Larch_C[[i]]<- grep("-",Larch[i,], value = TRUE, fixed = TRUE)
}


#subset to only rapidly expanding/ contracting fams (ones that contain "*")
Suillus_only_C_rapid <- list()
for (i in 1:length(Suillus_only_C)){
  Suillus_only_C_rapid[[i]]<- grep("*",Suillus_only_C[[i]], value = TRUE, fixed = TRUE)
}

Suillus_only_E_rapid <- list()
for (i in 1:length(Suillus_only_E)){
  Suillus_only_E_rapid[[i]]<- grep("*",Suillus_only_E[[i]], value = TRUE, fixed = TRUE)
}

Suillus_none_C_rapid <- list()
for (i in 1:length(Suillus_none_C)){
  Suillus_none_C_rapid[[i]]<- grep("*",Suillus_none_C[[i]], value = TRUE, fixed = TRUE)
}

Suillus_none_E_rapid <- list()
for (i in 1:length(Suillus_none_E)){
  Suillus_none_E_rapid[[i]]<- grep("*",Suillus_none_E[[i]], value = TRUE, fixed = TRUE)
}

Red_C_rapid <- list()
for (i in 1:length(Red_C)){
  Red_C_rapid[[i]]<- grep("*",Red_C[[i]], value = TRUE, fixed = TRUE)
}

Red_E_rapid <- list()
for (i in 1:length(Red_E)){
  Red_E_rapid[[i]]<- grep("*",Red_E[[i]], value = TRUE, fixed = TRUE)
}
White_C_rapid <- list()
for (i in 1:length(White_C)){
  White_C_rapid[[i]]<- grep("*",White_C[[i]], value = TRUE, fixed = TRUE)
}

White_E_rapid <- list()
for (i in 1:length(White_E)){
  White_E_rapid[[i]]<- grep("*",White_E[[i]], value = TRUE, fixed = TRUE)
}
Larch_C_rapid <- list()
for (i in 1:length(Larch_C)){
  Larch_C_rapid[[i]]<- grep("*",Larch_C[[i]], value = TRUE, fixed = TRUE)
}

Larch_E_rapid <- list()
for (i in 1:length(Larch_E)){
  Larch_E_rapid[[i]]<- grep("*",Larch_E[[i]], value = TRUE, fixed = TRUE)
}


#refomat as data frame
Suillus_only_E_df<- as.data.frame(unlist(Suillus_only_E_rapid))
Suillus_only_C_df<- as.data.frame(unlist(Suillus_only_C_rapid))
Suillus_none_E_df<- as.data.frame(unlist(Suillus_none_E_rapid))
Suillus_none_C_df<- as.data.frame(unlist(Suillus_none_C_rapid))

Red_C_df<- as.data.frame(unlist(Red_C_rapid))
Red_E_df<- as.data.frame(unlist(Red_E_rapid))
White_C_df<- as.data.frame(unlist(White_C_rapid))
White_E_df<- as.data.frame(unlist(White_E_rapid))
Larch_C_df<- as.data.frame(unlist(Larch_C_rapid))
Larch_E_df<- as.data.frame(unlist(Larch_E_rapid))

#clean for printing 
#remove numbers and brackets 
#, fixed = TRUE
Suillus_E_clean<- as.data.frame(lapply(Suillus_only_E_df, gsub, pattern = "\\[.*", replacement = "", perl=TRUE))
Suillus_C_clean<- as.data.frame(lapply(Suillus_only_C_df, gsub, pattern = "\\[.*", replacement = "", perl=TRUE))
Suillus_none_E_clean<- as.data.frame(lapply(Suillus_none_E_df, gsub, pattern = "\\[.*", replacement = "", perl=TRUE))
Suillus_none_C_clean<- as.data.frame(lapply(Suillus_none_C_df, gsub, pattern = "\\[.*", replacement = "", perl=TRUE))

Red_E_clean<- as.data.frame(lapply(Red_E_df, gsub, pattern = "\\[.*", replacement = "", perl=TRUE))
Red_C_clean<- as.data.frame(lapply(Red_C_df, gsub, pattern = "\\[.*", replacement = "", perl=TRUE))
White_E_clean<- as.data.frame(lapply(White_E_df, gsub, pattern = "\\[.*", replacement = "", perl=TRUE))
White_C_clean<- as.data.frame(lapply(White_C_df, gsub, pattern = "\\[.*", replacement = "", perl=TRUE))
Larch_E_clean<- as.data.frame(lapply(Larch_E_df, gsub, pattern = "\\[.*", replacement = "", perl=TRUE))
Larch_C_clean<- as.data.frame(lapply(Larch_C_df, gsub, pattern = "\\[.*", replacement = "", perl=TRUE))


#extract only the family number (line number in MCL file)
Suillus_E_cleaner<- as.data.frame(lapply(Suillus_E_clean, gsub, pattern = "MCLFAM", replacement = "", perl=TRUE))
Suillus_C_cleaner<- as.data.frame(lapply(Suillus_C_clean, gsub, pattern = "MCLFAM", replacement = "", perl=TRUE))
Suillus_none_E_cleaner<- as.data.frame(lapply(Suillus_none_E_clean, gsub, pattern = "MCLFAM", replacement = "", perl=TRUE))
Suillus_none_C_cleaner<- as.data.frame(lapply(Suillus_none_C_clean, gsub, pattern = "MCLFAM", replacement = "", perl=TRUE))

Red_E_cleaner<- as.data.frame(lapply(Red_E_clean, gsub, pattern = "MCLFAM", replacement = "", perl=TRUE))
Red_C_cleaner<- as.data.frame(lapply(Red_C_clean, gsub, pattern = "MCLFAM", replacement = "", perl=TRUE))
White_E_cleaner<- as.data.frame(lapply(White_E_clean, gsub, pattern = "MCLFAM", replacement = "", perl=TRUE))
White_C_cleaner<- as.data.frame(lapply(White_C_clean, gsub, pattern = "MCLFAM", replacement = "", perl=TRUE))
Larch_E_cleaner<- as.data.frame(lapply(Larch_E_clean, gsub, pattern = "MCLFAM", replacement = "", perl=TRUE))
Larch_C_cleaner<- as.data.frame(lapply(Larch_C_clean, gsub, pattern = "MCLFAM", replacement = "", perl=TRUE))

Suillus_E_cleaner2<- as.data.frame(lapply(Suillus_E_cleaner, gsub, pattern = "(^|[^0-9])0+", "\\1", replacement = "", perl=TRUE))
Suillus_C_cleaner2<- as.data.frame(lapply(Suillus_C_cleaner, gsub, pattern = "(^|[^0-9])0+", "\\1", replacement = "", perl=TRUE))
Suillus_none_E_cleaner2<- as.data.frame(lapply(Suillus_none_E_cleaner, gsub, pattern = "(^|[^0-9])0+", "\\1", replacement = "", perl=TRUE))
Suillus_none_C_cleaner2<- as.data.frame(lapply(Suillus_none_C_cleaner, gsub, pattern = "(^|[^0-9])0+", "\\1", replacement = "", perl=TRUE))

Red_E_cleaner2<- as.data.frame(lapply(Red_E_cleaner, gsub, pattern = "(^|[^0-9])0+", "\\1", replacement = "", perl=TRUE))
Red_C_cleaner2<- as.data.frame(lapply(Red_C_cleaner, gsub, pattern = "(^|[^0-9])0+", "\\1", replacement = "", perl=TRUE))
White_E_cleaner2<- as.data.frame(lapply(White_E_cleaner, gsub, pattern = "(^|[^0-9])0+", "\\1", replacement = "", perl=TRUE))
White_C_cleaner2<- as.data.frame(lapply(White_C_cleaner, gsub, pattern = "(^|[^0-9])0+", "\\1", replacement = "", perl=TRUE))
Larch_E_cleaner2<- as.data.frame(lapply(Larch_E_cleaner, gsub, pattern = "(^|[^0-9])0+", "\\1", replacement = "", perl=TRUE))
Larch_C_cleaner2<- as.data.frame(lapply(Larch_C_cleaner, gsub, pattern = "(^|[^0-9])0+", "\\1", replacement = "", perl=TRUE))


#get unique fams 
Suillus_E_cleaner3<- unique(Suillus_E_cleaner2)
Suillus_C_cleaner3<- unique(Suillus_C_cleaner2)
Suillus_none_E_cleaner3<- unique(Suillus_none_E_cleaner2)
Suillus_none_C_cleaner3<- unique(Suillus_none_C_cleaner2)

Red_E_cleaner3<- unique(Red_E_cleaner2)
Red_C_cleaner3<- unique(Red_C_cleaner2)
White_E_cleaner3<- unique(White_E_cleaner2)
White_C_cleaner3<- unique(White_C_cleaner2)
Larch_E_cleaner3<- unique(Larch_E_cleaner2)
Larch_C_cleaner3<- unique(Larch_C_cleaner2)


#call table 
Suillus_E_table<- as.data.frame(sort(table(Suillus_E_cleaner2), decreasing = TRUE))
Suillus_C_table<- as.data.frame(sort(table(Suillus_C_cleaner2), decreasing = TRUE))
Suillus_none_E_table<- as.data.frame(sort(table(Suillus_none_E_cleaner2), decreasing = TRUE))
Suillus_none_C_table<- as.data.frame(sort(table(Suillus_none_C_cleaner2), decreasing = TRUE))

Red_E_table<- as.data.frame(sort(table(Red_E_cleaner2), decreasing = TRUE))
Red_C_table<- as.data.frame(sort(table(Red_C_cleaner2), decreasing = TRUE))
White_E_table<- as.data.frame(sort(table(White_E_cleaner2), decreasing = TRUE))
White_C_table<- as.data.frame(sort(table(White_C_cleaner2), decreasing = TRUE))
Larch_E_table<- as.data.frame(sort(table(Larch_E_cleaner2), decreasing = TRUE))
Larch_C_table<- as.data.frame(sort(table(Larch_C_cleaner2), decreasing = TRUE))


#expansions in Suillus, not in others
other_fams_expanded<- Suillus_none_E_table$Suillus_none_E_cleaner2
Suillus_exclsuive_expansions<- Suillus_E_table[! Suillus_E_table$Suillus_E_cleaner2 %in% other_fams_expanded,]
View(Suillus_exclsuive_expansions)

#contractions in Suillus, not in others
other_fams_contracted<-Suillus_none_C_table$Suillus_none_C_cleaner2
Suillus_exclsuive_contractions<- Suillus_C_table[! Suillus_C_table$Suillus_C_cleaner2 %in% other_fams_contracted,]




#drop levels
White_C_table$White_C_cleaner2<- levels(droplevels(White_C_table$White_C_cleaner2))
Red_C_table$Red_C_cleaner2<- levels(droplevels(Red_C_table$Red_C_cleaner2))
Larch_C_table$Larch_C_cleaner2<- levels(droplevels(Larch_C_table$Larch_C_cleaner2))


#do this for R/W/L as well- are there rapid fam expansions / contractions that are specific to R/W/L?
#red
white_and_larch_contracted<- c(White_C_table$White_C_cleaner2, Larch_C_table$Larch_C_cleaner2)
white_and_larch_contracted<-unique(white_and_larch_contracted)
Red_exclsuive_contractions<-Red_C_table[! Red_C_table$Red_C_cleaner2 %in% white_and_larch_contracted,]

#white
red_and_larch_contracted<- c(Red_C_table$Red_C_cleaner2, Larch_C_table$Larch_C_cleaner2)
red_and_larch_contracted<-unique(red_and_larch_contracted)
White_exclsuive_contractions<-White_C_table[! White_C_table$White_C_cleaner2 %in% red_and_larch_contracted,]

#larch
red_and_white_contracted<- c(Red_C_table$Red_C_cleaner2, White_C_table$White_C_cleaner2)
red_and_white_contracted<-unique(red_and_white_contracted)
Larch_exclsuive_contractions<-Larch_C_table[! Larch_C_table$Larch_C_cleaner2 %in% red_and_white_contracted,]


#drop levels for expansions
White_E_table$White_E_cleaner2<- levels(droplevels(White_E_table$White_E_cleaner2))
Red_E_table$Red_E_cleaner2<- levels(droplevels(Red_E_table$Red_E_cleaner2))
Larch_E_table$Larch_E_cleaner2<- levels(droplevels(Larch_E_table$Larch_E_cleaner2))

#red
white_and_larch_expanded<- c(White_E_table$White_E_cleaner2, Larch_E_table$Larch_E_cleaner2)
white_and_larch_expanded<-unique(white_and_larch_expanded)
Red_exclsuive_expansions<-Red_E_table[! Red_E_table$Red_E_cleaner2 %in% white_and_larch_expanded,]

#white
red_and_larch_expanded<- c(Red_E_table$Red_E_cleaner2, Larch_E_table$Larch_E_cleaner2)
red_and_larch_expanded<-unique(red_and_larch_expanded)
White_exclsuive_expansions<-White_E_table[! White_E_table$White_E_cleaner2 %in% red_and_larch_expanded,]

#larch
red_and_white_expanded<- c(Red_E_table$Red_E_cleaner2, White_E_table$White_E_cleaner2)
red_and_white_expanded<-unique(red_and_white_expanded)
Larch_exclsuive_expansions<-Larch_E_table[! Larch_E_table$Larch_E_cleaner2 %in% red_and_white_expanded,]

#not that interesting - none are conserved across all host associations. 

#sort
Suillus_E_cleaner3$unlist.Suillus_only_E_rapid.<- sort(as.numeric(Suillus_E_cleaner3$unlist.Suillus_only_E_rapid.))
Suillus_C_cleaner3$unlist.Suillus_only_C_rapid.<- sort(as.numeric(Suillus_C_cleaner3$unlist.Suillus_only_C_rapid.))

#print these gene lists and extract in bash 
write.table(Suillus_E_cleaner3, "Suillus_E_line_numbers.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(Suillus_C_cleaner3, "Suillus_C_line_numbers.txt", sep="\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#get gene names for the families that are significantly expanded or contracted.
#MCL<-  read.table("All_MCL_Fams.I15.filtered.fixed.tab", header = FALSE, sep = "\t", check.names = FALSE, fill = TRUE)
MCL<-  read.table("dump.out.All.I15", header = FALSE, sep = "\t", check.names = FALSE, fill = TRUE)

#sort
Suillus_E_ordered<- sort(as.numeric(Suillus_E_cleaner3$unlist.Suillus_only_E_rapid.))
Suillus_C_ordered<- sort(as.numeric(Suillus_C_cleaner3$unlist.Suillus_only_C_rapid.))

expanded_genes_in_ea_fam<- MCL[Suillus_E_ordered, ]
row.names(expanded_genes_in_ea_fam)<- Suillus_E_ordered
contracted_genes_in_ea_fam<- MCL[Suillus_C_ordered, ]
row.names(contracted_genes_in_ea_fam)<- Suillus_C_ordered


#print off the gene files of interst for each family expansion / contraction 
write.table(expanded_genes_in_ea_fam, "expanded_genes_in_ea_fam.txt", sep="\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
write.table(contracted_genes_in_ea_fam, "contracted_genes_in_ea_fam.txt", sep="\t", row.names = TRUE, col.names = FALSE, quote = FALSE)
