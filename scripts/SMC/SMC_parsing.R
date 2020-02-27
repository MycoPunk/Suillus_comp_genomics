
#data was generated using fungal antiSMASH
#setwd("~/Desktop/Project_Suillus_comp_genomics/R")
setwd("~/Desktop/Project_Suillus_comp_genomics/Suillus_comp_genomics/data")

#inputs were genwerated from AS5 and parsed using the follosing script 
#fung<- data.frame(read.csv("<sp>_AS5.csv", header = TRUE))
#table(fung$class)
#table<- table(fung$class)
#sum(table)

#load packages 
library(tidyr)
library(reshape2)
library(dplyr)
library(car)
library(ggplot2)
install.packages("cfcdae_0.8-4.tgz",repos=NULL)
library("cfcdae")
library("lme4")

#read in the input files
SMCs_DF<- read.csv("AS5_outputs.csv")

#shrink DF to just SMCs
SMC_DF_small<- SMCs_DF[,c(2, 5:length(SMCs_DF))]

#check for zeros and remove if found
colSums(SMC_DF_small[2:length(SMC_DF_small)])
#there are no columns that sum to zero 
length(colnames(SMC_DF_small))
#there are 13 categories of SMCs

##seperate the two groups
#shrink DF to the important stuff
SMCs_DF_Suillus<- SMCs_DF[SMCs_DF$treatment == "S",]
SMCs_DF_Suillus2<- SMCs_DF_Suillus[,5:length(SMCs_DF_Suillus)]
#get means for each SMC type
SMCs_DF_Suillus3<- data.frame(round(colMeans(SMCs_DF_Suillus2),0))
#reformat into df with the correct row and col names
names<- rownames(SMCs_DF_Suillus3)
SMCs_DF_Suillus4<- cbind(SMCs_DF_Suillus3, names, rep(("a_Suillus"), length(SMCs_DF_Suillus3)))
colnames(SMCs_DF_Suillus4)<- c("counts", "names", "group")

SMCs_DF_Other<- SMCs_DF[SMCs_DF$treatment == "O",]
#shrink DF to the important stuff
SMCs_DF_Other2<- SMCs_DF_Other[,5:length(SMCs_DF_Other)]
#get means for each SMC type
SMCs_DF_Other3<- data.frame(round(colMeans(SMCs_DF_Other2),0))
#reformat into df with the correct row and col names
names<- rownames(SMCs_DF_Other3)
SMCs_DF_Other4<- cbind(SMCs_DF_Other3, names, rep(("b_Other"), length(SMCs_DF_Other3)))
colnames(SMCs_DF_Other4)<- c("counts", "names", "group")

#bind them 
SMCs_DF<- rbind(SMCs_DF_Suillus4, SMCs_DF_Other4)
#make repeating rows by SMC counts
counts_gather2 <- as.data.frame(lapply(SMCs_DF, rep, SMCs_DF$counts))
#wants to inherit zeros in the following- get rid of them by switching types (I dono)
counts_gather3<- as.matrix(counts_gather2)
counts_gather4<- as.data.frame(counts_gather3)

#call table
counts_gather_table<- table(counts_gather4$group, counts_gather4$names)
groups<- c("Suillus", "Other")

#make spine plot
par(las = 1)
par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))
spineplot(counts_gather_table, main = " ", 
          col=c("#5676A1", "#405952","#969574", "#FDD191", "#885053", "#442B47"), 
          border = NA, 
          #cex.axis = 0.66, 
          #ylab = "", #xlab = "", 
          yaxlabels = NA, 
          xaxlabels = groups)


#legend for spineplot
legend("topleft", legend = colnames(counts_gather_table) , 
       fill=c("#5676A1", "#405952","#969574", "#FDD191", "#885053", "#442B47"),
       ncol = 1, 
       cex = .4, 
       lwd = 0, 
       box.lwd = 0, 
       box.lty =0, 
       box.col =0, 
       xjust = 1)


#get n and N
counts_gather_table
N<- rowSums(counts_gather_table)
N

##make grouped boxplots
#seperate Suillus and Other
counts_gather2_Suillus<- counts_gather2[counts_gather2$group == "a_Suillus",]
counts_gather2_Other<- counts_gather2[counts_gather2$group == "b_Other",]

#put data into long form
#add repeat col (ours is only 1)
SMC_DF_small$id<- seq(1:nrow(SMC_DF_small))
#organize
SMC_DF_small_long <- SMC_DF_small %>% group_by(treatment) %>%
  gather(data = SMC_DF_small, id, T1PKS, fungal.RiPP, terpene.T1PKS, Terpene, indole, NRPS.like.siderophore, betalactone, NRPS, sideraphore, NRPS.like.Terpene, NRPS.like, NRPS.T1PKS, NRPS.like.T1PKS)

#make col names match
colnames(SMC_DF_small_long)<- c("Group", "SMC", "Count")

#subset to the primary groups (> 1 average occurance for that group)
names_primary_groups<- colnames(counts_gather_table)
primary_only<- SMC_DF_small_long[SMC_DF_small_long$SMC %in% names_primary_groups,]

#make grouped box plot
p<- ggplot(primary_only, aes(x=SMC, y=Count, fill=Group)) + 
  geom_boxplot()
p+scale_fill_manual(values=c("#B1913A", "#273666")) + theme_minimal()

#add points
p + geom_point(position=position_jitterdodge(jitter.width = .13),alpha=0.3) +
  theme_bw(base_size = 3) + scale_fill_manual(values=c("#B1913A", "#273666")) + theme_minimal()

#horizontal
p + geom_point(position=position_jitterdodge(jitter.width = .13),alpha=0.3) +
  theme_bw(base_size = 3) + scale_fill_manual(values=c("#B1913A", "#273666")) + theme_minimal() + coord_flip()



#####stats start here 
#####ANOVA 
#type 2 ANOVA

#make model
model.for.t2<-lm(primary_only$Count ~ primary_only$SMC * primary_only$Group)
model.for.t3<- Anova(model.for.t2, type=2)
model.for.t3

#run TukeyHSD
#define interaction term
primary_only$Group_SMC_interaction <- as.factor(interaction(primary_only$SMC, primary_only$Group))

#re-run model with interaction term
model.for.t4<- lm(primary_only$Count ~ primary_only$Group_SMC_interaction)

#take a look - lines in common mean not significantly different 
sidelines(pairwise(model.for.t4, primary_only$Group_SMC_interaction,confidence = 0.95, type = "hsd"))
lines(pairwise(model.for.t4, primary_only$Group_SMC_interaction,confidence = 0.95, type = "hsd"))

#full result with hsd -- note you can't get p-vals with hsd though with cfcdae. 
pairwise(model.for.t4, primary_only$Group_SMC_interaction,confidence = 0.95, type = "hsd")

#run model for multiple t tests so that you can get p-vals
pairs_out<- pairwise(model.for.t4, primary_only$Group_SMC_interaction,confidence = 0.95, type = "regwr")

lines(pairs_out)

#get p-value using mulitple t tests and holm adjustment for multiple comparisons. 
p_vals<-with(primary_only, pairwise.t.test(primary_only$Count, primary_only$Group_SMC_interaction,
                                     p.adjust.method="holm"))


#make into dataframe to view
easy_view<- as.data.frame(p_vals[3])
#T1PKS S-O = 1.000000e+00
#terpene S-O = 6.939993e-19
#indole S-O = 1.000000e+00
#NRPS S-O = 1.000000e+00
#NRPS.like.Terpene S-O = 1.000000e+00
#NRPS.like S-O = 2.516576e-06


###### Within Suillus by host comparison ######

#parse dataframe by host association
#red
red_SMC_df<- SMCs_DF_Suillus[SMCs_DF_Suillus$host == "R",]
#white
white_SMC_df<- SMCs_DF_Suillus[SMCs_DF_Suillus$host == "W",]
#larch
larch_SMC_df<- SMCs_DF_Suillus[SMCs_DF_Suillus$host == "L",]


#format dataframs to relevent info
#red
red_SMC_df2<- red_SMC_df[,5:length(red_SMC_df)]
red_SMC_df3<- data.frame(round(colMeans(red_SMC_df2),0))
names<- rownames(red_SMC_df3)
red_SMC_df4<- cbind(red_SMC_df3, names, rep(("a_Red"), length(red_SMC_df3)))
colnames(red_SMC_df4)<- c("counts", "names", "group")

#white
white_SMC_df2<- white_SMC_df[,5:length(white_SMC_df)]
white_SMC_df3<- data.frame(round(colMeans(white_SMC_df2),0))
names<- rownames(white_SMC_df3)
white_SMC_df4<- cbind(white_SMC_df3, names, rep(("b_White"), length(white_SMC_df3)))
colnames(white_SMC_df4)<- c("counts", "names", "group")

#larch
larch_SMC_df2<- larch_SMC_df[,5:length(larch_SMC_df)]
larch_SMC_df3<- data.frame(round(colMeans(larch_SMC_df2),0))
names<- rownames(larch_SMC_df3)
larch_SMC_df4<- cbind(larch_SMC_df3, names, rep(("c_Larch"), length(larch_SMC_df3)))
colnames(larch_SMC_df4)<- c("counts", "names", "group")

#bind them 
SMCs_Suillus_DF<- rbind(red_SMC_df4, white_SMC_df4, larch_SMC_df4)

#shrink to ony the primary SMC groups
primary_suillus_only<- SMCs_Suillus_DF[SMCs_Suillus_DF$names %in% names_primary_groups,]

#make repeating rows by SMC counts
counts_gather2 <- as.data.frame(lapply(primary_suillus_only, rep, primary_suillus_only$counts))
#wants to inherit zeros in the following- get rid of them by switching types (I dono)
counts_gather3<- as.matrix(counts_gather2)
counts_gather4<- as.data.frame(counts_gather3)

#call table
counts_gather_table<- table(counts_gather4$group, counts_gather4$names)

groups<- c("Red", "White", "Larch")
#spine plot - note - I removed hte light green color for "indole" as there are not indols in the Suillus group 
par(las = 1)
par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))
spineplot(counts_gather_table, main = " ", 
          col=c("#405952","#969574", "#FDD191", "#885053", "#442B47"), 
          border = NA, 
          #cex.axis = 0.66, 
          #ylab = "", #xlab = "", 
          yaxlabels = NA, 
          xaxlabels = groups)


#legend for spineplot
legend("topleft", legend = colnames(counts_gather_table) , 
       fill=c("#405952","#969574", "#FDD191", "#885053", "#442B47"),
       ncol = 1, 
       cex = .4, 
       lwd = 0, 
       box.lwd = 0, 
       box.lty =0, 
       box.col =0, 
       xjust = 1)


#get n and N
counts_gather_table
N<- rowSums(counts_gather_table)
N


####make grouped box plot
#first put into long form
SMCs_DF_Suillus2<- SMCs_DF_Suillus[c(SMCs_DF_Suillus$host == "W" | SMCs_DF_Suillus$host == "R" | SMCs_DF_Suillus$host == "L"),]

#shrink and remove cols == zero
#only SMCs in primary list:
SMC_DF_small_no_zeros_Suillus<- SMCs_DF_Suillus2[,colnames(SMCs_DF_Suillus2) %in% names_primary_groups]
row.names(SMC_DF_small_no_zeros_Suillus)<- SMCs_DF_Suillus2$species
#cols not equal to zero: (indole)
SMC_DF_small_no_zeros_Suillus2<- SMC_DF_small_no_zeros_Suillus[,colSums(SMC_DF_small_no_zeros_Suillus) > 0]
#shrink to only relevent cols

#add back in host
SMC_DF_small_no_zeros_Suillus2$host<- SMCs_DF_Suillus2$host
#add repeat col (ours is only 1)
SMC_DF_small_no_zeros_Suillus2$id<- seq(1:nrow(SMC_DF_small_no_zeros_Suillus))
SMC_DF_small_no_zeros_Suillus2 <- SMC_DF_small_no_zeros_Suillus2 %>% group_by(host) %>%
  gather(data = SMC_DF_small_no_zeros_Suillus2, id, NRPS, NRPS.like, NRPS.like.Terpene, T1PKS, Terpene)

#make col names match
colnames(SMC_DF_small_no_zeros_Suillus2)<- c("Group", "SMC", "Count")


#set order of boxes
SMC_DF_small_no_zeros_Suillus2$Group <- factor(SMC_DF_small_no_zeros_Suillus2$Group,
                       levels = c('L','W', 'R'),ordered = TRUE)

### make grouped box plot
p<- ggplot(SMC_DF_small_no_zeros_Suillus2, aes(x=SMC, y=Count, fill=Group)) + 
  geom_boxplot()
p+scale_fill_manual(values=c("#B1913A", "#273666", "#FCD090")) + theme_minimal()

#horizontal
p + geom_point(position=position_jitterdodge(jitter.width = .13),alpha=0.3) +
  theme_bw(base_size = 3) + scale_fill_manual(values=c("#B1913A", "#273666", "#FCD090")) + theme_minimal() + coord_flip()


####### Stats #######
#type 2 ANOVA
#make model
model.for.t2_Suillus_only<-lm(SMC_DF_small_no_zeros_Suillus2$Count ~ SMC_DF_small_no_zeros_Suillus2$SMC * SMC_DF_small_no_zeros_Suillus2$Group)

model.for.t3<- Anova(model.for.t2_Suillus_only, type=2)
model.for.t3
#no significant difference between groups. 
