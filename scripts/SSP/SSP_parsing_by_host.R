
#load libraries
library("seqinr")
library("data.table")
library("stringr")
library("car")
library("FSA")
library("rcompanion")
library("coin")
library("multcompView")
library("multcomp")

#to load coin
#if(!require(coin)){install.packages("coin")}
#if(!require(FSA)){install.packages("FSA")}
#if(!require(rcompanion)){install.packages("rcompanion")}
#if(!require(multcompView)){install.packages("multcompView")}

options(stringsAsFactors = FALSE)
setwd("~/Desktop/Project_Suillus_comp_genomics/R")

#read it back in 
#SSPs_coded_within_Suillus_temp<- read.csv("SSPs_coded_within_Suillus.csv")
SSPs_coded_within_Suillus<- read.csv("SSP_and_Effector_totals.csv")


#split the two data categories for the three groups and get stats for each 
Red_df<- SSPs_coded_within_Suillus[SSPs_coded_within_Suillus$group=="a_Red",]
White_df<- SSPs_coded_within_Suillus[SSPs_coded_within_Suillus$group=="b_White",]
Larch_df<- SSPs_coded_within_Suillus[SSPs_coded_within_Suillus$group=="c_Larch",]

#check that the seperation worked
nrow(Red_df)
nrow(White_df)
nrow(Larch_df)
#yep

#means
Red_genome_size_mean<- mean(Red_df$genome_size)
White_genome_size_mean<- mean(White_df$genome_size)
Larch_genome_size_mean<- mean(Larch_df$genome_size)

#SD and SE
std <- function(x) sd(x)/sqrt(length(x))

Red_genome_size_SD<- sd(Red_df$genome_size)
Red_genome_size_SE<- std(Red_df$genome_size)  

White_genome_size_SD<- sd(White_df$genome_size)
White_genome_size_SE<- std(White_df$genome_size)  

Larch_genome_size_SD<- sd(Larch_df$genome_size)
Larch_genome_size_SE<- std(Larch_df$genome_size)  

all_df<- rbind(Red_df, White_df, Larch_df)

#build model
genome_size_model <- lm(all_df$genome_size ~ all_df$group)

#test normality and var
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(genome_size_model) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
boxCox(genome_size_model, lambda = seq(-5, 1, 1/10)) #what does boxcox suggest? log transformation works. 

#log transform to improve normality?
all_df$genome_size_log<- log(all_df$genome_size)

#see if the ransformation helped
genome_size_model2 <- lm(all_df$genome_size_log ~ all_df$group)
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(genome_size_model2) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
#not a lot, try -2 ==(1/x^2) - inverse square transformation 

all_df$genome_size_inv_suqare<- all_df$genome_size^(-2)

genome_size_model3 <- lm(all_df$genome_size_inv_suqare ~ all_df$group)
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(genome_size_model3) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
#that looks a little better 

#Run anova
summary(aov(genome_size_model3))
#not significant

####get p val with a randomization test
#run permutation test with coin for three factors


#make df
r_GS<- data.frame(Red_df$genome_size, rep("Red"))
colnames(r_GS)<- c("GS", "group")
w_GS<- data.frame(White_df$genome_size, rep("White"))
colnames(w_GS)<- c("GS", "group")
l_GS<- data.frame(Larch_df$genome_size, rep("Larch"))
colnames(l_GS)<- c("GS", "group")
#bind them
GS_df<- rbind(r_GS, w_GS, l_GS)
#set factors
#set order of boxes
GS_df$group <- factor(GS_df$group,levels = c('Red','White', 'Larch'),ordered = TRUE)

#run independence test (this is like an ANOVA)
independence_test(GS ~ group, 
                  data = GS_df)



####proteins
#means
Red_prot_mean<- mean(Red_df$n_proteins)
White_prot_mean<- mean(White_df$n_proteins)
Larch_prot_mean<- mean(Larch_df$n_proteins)

#SD and SE
Red_prot_SD<- sd(Red_df$n_proteins)
Red_prot_SE<- std(Red_df$n_proteins)  

White_prot_SD<- sd(White_df$n_proteins)
White_prot_SE<- std(White_df$n_proteins)  

Larch_prot_SD<- sd(Larch_df$n_proteins)
Larch_prot_SE<- std(Larch_df$n_proteins)  

#build model
prot_model <- lm(all_df$n_proteins ~ all_df$group)

#vis normality and variance
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(prot_model) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
boxCox(prot_model) #what does boxcox suggest? log transformation works. 
powerTransform(prot_model) #how about powerTransform? #log should be fine for this.

#log transform to improve normality?
all_df$n_proteins_log<- log(all_df$n_proteins)

#see if the transformation helped
prot_model2 <- lm(all_df$n_proteins_log ~ all_df$group)
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(prot_model2) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
#well, that didn't help...

#try the boc recomended transformation
prot_model3 <- lm(all_df$n_proteins_cat^-.5 ~ all_df$group)
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(prot_model3) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
#worse! go back. do not pass go. do not colelct $100. Use original data. 

#run anova
summary(aov(prot_model))
#not significant

#make df
r_GS<- data.frame(Red_df$n_proteins, rep("Red"))
colnames(r_GS)<- c("GS", "group")
w_GS<- data.frame(White_df$n_proteins, rep("White"))
colnames(w_GS)<- c("GS", "group")
l_GS<- data.frame(Larch_df$n_proteins, rep("Larch"))
colnames(l_GS)<- c("GS", "group")
#bind them
GS_df<- rbind(r_GS, w_GS, l_GS)
#set factors
#set order of boxes
GS_df$group <- factor(GS_df$group,levels = c('Red','White', 'Larch'),ordered = TRUE)

#run independence test (this is like an ANOVA)
independence_test(GS ~ group, 
                  data = GS_df)




####SSPs
#means
Red_ssp_mean<- mean(Red_df$n_SSPs)
White_ssp_mean<- mean(White_df$n_SSPs)
Larch_ssp_mean<- mean(Larch_df$n_SSPs)

#SD and SE
Red_ssp_SD<- sd(Red_df$n_SSPs)
Red_ssp_SE<- std(Red_df$n_SSPs)  

White_ssp_SD<- sd(White_df$n_SSPs)
White_ssp_SE<- std(White_df$n_SSPs)  

Larch_ssp_SD<- sd(Larch_df$n_SSPs)
Larch_ssp_SE<- std(Larch_df$n_SSPs)  

#build model
ssp_model <- lm(all_df$n_SSPs ~ all_df$group)

#vis normality and variance
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(ssp_model) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
boxCox(ssp_model, lambda = seq(-6,2,2)) #what does boxcox suggest?  
powerTransform(ssp_model) #how about powerTransform? 
# log to ^-2 

all_df$n_SSPs_log<- log(all_df$n_SSPs)

ssp_model2 <- lm(all_df$n_SSPs_log ~ all_df$group)

#see if the transformation helped
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(ssp_model2) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
#looks kinda better

#run anova
summary(aov(ssp_model2))
#not significant


#make df
r_GS<- data.frame(Red_df$n_SSPs, rep("Red"))
colnames(r_GS)<- c("GS", "group")
w_GS<- data.frame(White_df$n_SSPs, rep("White"))
colnames(w_GS)<- c("GS", "group")
l_GS<- data.frame(Larch_df$n_SSPs, rep("Larch"))
colnames(l_GS)<- c("GS", "group")
#bind them
GS_df<- rbind(r_GS, w_GS, l_GS)
#set factors
#set order of boxes
GS_df$group <- factor(GS_df$group,levels = c('Red','White', 'Larch'),ordered = TRUE)

#run independence test (this is like an ANOVA)
independence_test(GS ~ group, 
                  data = GS_df)



####SSSPs
#means
Red_ssp_mean<- mean(Red_df$nSSSPs)
White_ssp_mean<- mean(White_df$nSSSPs)
Larch_ssp_mean<- mean(Larch_df$nSSSPs)

#SD and SE
Red_sssp_SD<- sd(Red_df$nSSSPs)
Red_sssp_SE<- std(Red_df$nSSSPs)  

White_sssp_SD<- sd(White_df$nSSSPs)
White_sssp_SE<- std(White_df$nSSSPs)  

Larch_sssp_SD<- sd(Larch_df$nSSSPs)
Larch_sssp_SE<- std(Larch_df$nSSSPs)  

#build model
sssp_model <- lm(all_df$nSSSPs ~ all_df$group)

#vis normality and variance
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(sssp_model) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
boxCox(sssp_model, lambda = seq(-6,2,2)) #what does boxcox suggest?  
powerTransform(sssp_model) #how about powerTransform? don't bother - it looks pretty good. 


#run anova
summary(aov(sssp_model))
#significant!

#what's different?
sssp_aov<- aov(sssp_model)
TukeyHSD(sssp_aov)
#red and larch are sig. different 

#randomization test make df
r_GS<- data.frame(Red_df$nSSSPs, rep("Red"))
colnames(r_GS)<- c("GS", "group")
w_GS<- data.frame(White_df$nSSSPs, rep("White"))
colnames(w_GS)<- c("GS", "group")
l_GS<- data.frame(Larch_df$nSSSPs, rep("Larch"))
colnames(l_GS)<- c("GS", "group")
#bind them
GS_df<- rbind(r_GS, w_GS, l_GS)
#set factors
#set order of boxes
GS_df$group <- factor(GS_df$group,levels = c('Red','White', 'Larch'),ordered = TRUE)

#run independence test (this is like an ANOVA)
independence_test(GS ~ group, 
                  data = GS_df)
#significant!

#run tukey equivelent to see whats significant
#take a look at the medians 
boxplot(GS ~ group,
        data = GS_df)

#order them by median (highest to lowest)
GS_df$group = factor(GS_df$group , 
                     levels = c("Larch", "White", "Red"))


headtail(GS_df)

### Pairwise tests


PT = pairwisePermutationTest(GS ~ group,
                             data = GS_df,
                             method="fdr")

PT


####Effectors
#means
Red_effector_mean<- mean(Red_df$n_Effectors)
White_effector_mean<- mean(White_df$n_Effectors)
Larch_effector_mean<- mean(Larch_df$n_Effectors)

#SD and SE
Red_effector_SD<- sd(Red_df$n_Effectors)
Red_effector_SE<- std(Red_df$n_Effectors)  

White_effector_SD<- sd(White_df$n_Effectors)
White_effector_SE<- std(White_df$n_Effectors)  

Larch_effector_SD<- sd(Larch_df$n_Effectors)
Larch_effector_SE<- std(Larch_df$n_Effectors)  

#build model
effector_model <- lm(all_df$n_Effectors ~ all_df$group)

#vis normality and variance
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(effector_model) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
boxCox(effector_model, lambda = seq(-6,2,2)) #what does boxcox suggest?  
powerTransform(effector_model) #how about powerTransform? zero is in there, try log transformation

#transform
all_df$n_Effectors_log<- log(all_df$n_Effectors)

#make new model
effector_model2 <- lm(all_df$n_Effectors_log ~ all_df$group)

#see if the transformation helped
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(effector_model2) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
#looks abotu the same, just skip the transformation

#run anova
summary(aov(effector_model))
#not significant

#randomization test make df
r_GS<- data.frame(Red_df$n_Effectors, rep("Red"))
colnames(r_GS)<- c("GS", "group")
w_GS<- data.frame(White_df$n_Effectors, rep("White"))
colnames(w_GS)<- c("GS", "group")
l_GS<- data.frame(Larch_df$n_Effectors, rep("Larch"))
colnames(l_GS)<- c("GS", "group")
#bind them
GS_df<- rbind(r_GS, w_GS, l_GS)
#set factors
#set order of boxes
GS_df$group <- factor(GS_df$group,levels = c('Red','White', 'Larch'),ordered = TRUE)

#run independence test (this is like an ANOVA)
independence_test(GS ~ group, 
                  data = GS_df)

#not significant





##### percentages analysis #####


#%SSPs out of all prot.
#means
Red_SSP_per_all_prot_mean<- mean(Red_df$percent_SSPs_out_of_total_genes)
White_SSP_per_all_prot_mean<- mean(White_df$percent_SSPs_out_of_total_genes)
Larch_SSP_per_all_prot_mean<- mean(Larch_df$percent_SSPs_out_of_total_genes)

#SD and SE
Red_SSP_per_all_prot_SD<- sd(Red_df$percent_SSPs_out_of_total_genes)
Red_SSP_per_all_prot_SE<- std(Red_df$percent_SSPs_out_of_total_genes)  

White_SSP_per_all_prot_SD<- sd(White_df$percent_SSPs_out_of_total_genes)
White_SSP_per_all_prot_SE<- std(White_df$percent_SSPs_out_of_total_genes)  

Larch_SSP_per_all_prot_SD<- sd(Larch_df$percent_SSPs_out_of_total_genes)
Larch_SSP_per_all_prot_SE<- std(Larch_df$percent_SSPs_out_of_total_genes) 

#build model
SSP_per_all_prot_model <- lm(all_df$percent_SSPs_out_of_total_genes ~ all_df$group)

#vis normality and variance
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(SSP_per_all_prot_model) #data is non normal, variance looks ok-ish
par(mfrow = c(1,1)) #back to one image per page
boxCox(SSP_per_all_prot_model, lambda = seq(-6,2,2)) #what does boxcox suggest?  
powerTransform(SSP_per_all_prot_model) #how about powerTransform? # ^-2 is what's recommended, but zero is in there, try log.

#transform
all_df$percent_SSPs_out_of_total_genes_log<- log(all_df$percent_SSPs_out_of_total_genes)

#make new model
SSP_per_all_prot_model2 <- lm(all_df$percent_SSPs_out_of_total_genes_log ~ all_df$group)

#see if the transformation helped
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(SSP_per_all_prot_model2) #data looks better
par(mfrow = c(1,1)) #back to one image per page
#looks about the same, just use the original


#run anova
summary(aov(SSP_per_all_prot_model))
#not significant

#randomization test make df
r_GS<- data.frame(Red_df$percent_SSPs_out_of_total_genes, rep("Red"))
colnames(r_GS)<- c("GS", "group")
w_GS<- data.frame(White_df$percent_SSPs_out_of_total_genes, rep("White"))
colnames(w_GS)<- c("GS", "group")
l_GS<- data.frame(Larch_df$percent_SSPs_out_of_total_genes, rep("Larch"))
colnames(l_GS)<- c("GS", "group")
#bind them
GS_df<- rbind(r_GS, w_GS, l_GS)
#set factors
#set order of boxes
GS_df$group <- factor(GS_df$group,levels = c('Red','White', 'Larch'),ordered = TRUE)

#run independence test (this is like an ANOVA)
independence_test(GS ~ group, 
                  data = GS_df)
#marganaly significant 

#whats different?
#run tukey equivelent to see whats significant
#take a look at the medians 
boxplot(GS ~ group,
        data = GS_df)

#order them by median (highest to lowest)
GS_df$group = factor(GS_df$group , 
                     levels = c("Larch", "White", "Red"))


headtail(GS_df)

### Pairwise tests
PT = pairwisePermutationTest(GS ~ group,
                             data = GS_df,
                             method="fdr")

PT

cldList(comparison = PT$Comparison,
        p.value    = PT$p.adjust,
        threshold  = 0.05)
#not significant after fdr correction


#%SSSPs out of SSPs.
#means
Red_SSSP_per_SSP_mean<- mean(Red_df$percent_SSSPs_out_of_SSPs)
White_SSSP_per_SSP_mean<- mean(White_df$percent_SSSPs_out_of_SSPs)
Larch_SSSP_per_SSP_mean<- mean(Larch_df$percent_SSSPs_out_of_SSPs)

#SD and SE
Red_SSSP_per_SSP_SD<- sd(Red_df$percent_SSSPs_out_of_SSPs)
Red_SSSP_per_SSP_SE<- std(Red_df$percent_SSSPs_out_of_SSPs)  

White_SSSP_per_SSP_SD<- sd(White_df$percent_SSSPs_out_of_SSPs)
White_SSSP_per_SSP_SE<- std(White_df$percent_SSSPs_out_of_SSPs)  

Larch_SSSP_per_SSP_SD<- sd(Larch_df$percent_SSSPs_out_of_SSPs)
Larch_SSSP_per_SSP_SE<- std(Larch_df$percent_SSSPs_out_of_SSPs) 

#build model
SSSP_per_SSP_model <- lm(all_df$percent_SSSPs_out_of_SSPs ~ all_df$group)

#vis normality and variance
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(SSSP_per_SSP_model) #data is non normal, variance looks pretty good
par(mfrow = c(1,1)) #back to one image per page
boxCox(SSSP_per_SSP_model) #what does boxcox suggest?  
powerTransform(SSSP_per_SSP_model) #recommended not to transform

#run anova
summary(aov(SSSP_per_SSP_model))
#significant! 

#what's different?
sssp_aov2<- aov(SSSP_per_SSP_model)
TukeyHSD(sssp_aov2)
#red and larch are sig. different 

#randomization test make df
r_GS<- data.frame(Red_df$percent_SSSPs_out_of_SSPs, rep("Red"))
colnames(r_GS)<- c("GS", "group")
w_GS<- data.frame(White_df$percent_SSSPs_out_of_SSPs, rep("White"))
colnames(w_GS)<- c("GS", "group")
l_GS<- data.frame(Larch_df$percent_SSSPs_out_of_SSPs, rep("Larch"))
colnames(l_GS)<- c("GS", "group")
#bind them
GS_df<- rbind(r_GS, w_GS, l_GS)
#set factors
#set order of boxes
GS_df$group <- factor(GS_df$group,levels = c('Red','White', 'Larch'),ordered = TRUE)

#run independence test (this is like an ANOVA)
independence_test(GS ~ group, 
                  data = GS_df)

#significant!


#run tukey equivelent to see whats significant
#take a look at the medians 
boxplot(GS ~ group,
        data = GS_df)

#order them by median (highest to lowest)
GS_df$group = factor(GS_df$group , 
                     levels = c("Larch", "White", "Red"))


headtail(GS_df)

### Pairwise tests


PT = pairwisePermutationTest(GS ~ group,
                             data = GS_df,
                             method="fdr")

PT




#Effectors out of SSPs
#means
Red_effectors_out_of_SSP_mean<- mean(Red_df$percent_effectors_out_of_SSPs)
White_effectors_out_of_SSP_mean<- mean(White_df$percent_effectors_out_of_SSPs)
Larch_effectors_out_of_SSP_mean<- mean(Larch_df$percent_effectors_out_of_SSPs)

#SD and SE
Red_effectors_out_of_SSP_SD<- sd(Red_df$percent_effectors_out_of_SSPs)
Red_effectors_out_of_SSP_SE<- std(Red_df$percent_effectors_out_of_SSPs)  

White_effectors_out_of_SSP_SD<- sd(White_df$percent_effectors_out_of_SSPs)
White_effectors_out_of_SSP_SE<- std(White_df$percent_effectors_out_of_SSPs)  

Larch_effectors_out_of_SSP_SD<- sd(Larch_df$percent_effectors_out_of_SSPs)
Larch_effectors_out_of_SSP_SE<- std(Larch_df$percent_effectors_out_of_SSPs) 

#build model
effectors_out_of_SSP_model <- lm(all_df$percent_effectors_out_of_SSPs ~ all_df$group)

#vis normality and variance
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(effectors_out_of_SSP_model) #data is non normal, with non-constant variance 
par(mfrow = c(1,1)) #back to one image per page
boxCox(effectors_out_of_SSP_model) #what does boxcox suggest?  
powerTransform(effectors_out_of_SSP_model) #how about powerTransform? # ^-3 is what's recommended, but -2 is in range. 

#transform
all_df$percent_effectors_out_of_SSPs_inv_square<- all_df$percent_effectors_out_of_SSPs^(-2)

#make new model
effectors_out_of_SSP_model2 <- lm(all_df$percent_effectors_out_of_SSPs_inv_square ~ all_df$group)

#see if the transformation helped
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(effectors_out_of_SSP_model2) #data looks better
par(mfrow = c(1,1)) #back to one image per page
#looks better enough 

#run anova
summary(aov(effectors_out_of_SSP_model2))
#not significant

#randomization test make df
r_GS<- data.frame(Red_df$percent_effectors_out_of_SSPs, rep("Red"))
colnames(r_GS)<- c("GS", "group")
w_GS<- data.frame(White_df$percent_effectors_out_of_SSPs, rep("White"))
colnames(w_GS)<- c("GS", "group")
l_GS<- data.frame(Larch_df$percent_effectors_out_of_SSPs, rep("Larch"))
colnames(l_GS)<- c("GS", "group")
#bind them
GS_df<- rbind(r_GS, w_GS, l_GS)
#set factors
#set order of boxes
GS_df$group <- factor(GS_df$group,levels = c('Red','White', 'Larch'),ordered = TRUE)

#run independence test (this is like an ANOVA)
independence_test(GS ~ group, 
                  data = GS_df)


######## plots start here ########
#plot of n_SSPs

#means from transformaton 
Red_df_new<- all_df[all_df$group=="a_Red",]
White_df_new<- all_df[all_df$group=="b_White",]
Larch_df_new<- all_df[all_df$group=="c_Larch",]

M_red<- mean(Red_df_new$n_SSPs)
M_white<- mean(White_df_new$n_SSPs)
M_larch<- mean(Larch_df_new$n_SSPs)

par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))
stripchart(all_df$n_SSPs ~ all_df$group,
           vertical = TRUE,
           method = "jitter", jitter = 0.2, 
           pch = 16, 
           col = c("#5676A166","#B0913666", "#B3675766"),
           bg = rep(c("#5676A166","#B0913666", "#B3675766")),
           cex.axis = 0.7,
           ylim=c(0,550), 
           ylab = expression(paste("SSPs")), 
           axes = FALSE, 
           cex = 1.3)
box()
axis(2)
mtext(text = c("Red", "White", "Larch"),side=1,at=c(1,2,3),line = 1, font = 1)
segments(x0 = .7, y0 =  M_red, x1 = 1.3, y1= M_red, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  M_white, x1 = 2.3, y1= M_white, lwd = 2, col = "black" )
segments(x0 = 2.7, y0 =  M_larch, x1 = 3.3, y1= M_larch, lwd = 2, col = "black" )
mtext(c("a", "a", "a"),side=1,at=c(1,2,3),line = -7, font = 1)


###plot of n SSSPs

#means from transformaton 
M_red<- mean(Red_df_new$nSSSPs)
M_white<- mean(White_df_new$nSSSPs)
M_larch<- mean(Larch_df_new$nSSSPs)

par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))
stripchart(all_df$nSSSPs ~ all_df$group,
           vertical = TRUE,
           method = "jitter", jitter = 0.2, 
           pch = 16, 
           col = c("#5676A166","#B0913666", "#B3675766"),
           bg = rep(c("#5676A166","#B0913666", "#B3675766")),
           cex.axis = 0.7,
           ylim=c(0,220), 
           ylab = expression(paste("SSSPs")), 
           axes = FALSE, 
           cex = 1.3)
box()
axis(2) 
mtext(text = c("Red", "White", "Larch"),side=1,at=c(1,2,3),line = 1, font = 1)
segments(x0 = .7, y0 =  M_red, x1 = 1.3, y1= M_red, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  M_white, x1 = 2.3, y1= M_white, lwd = 2, col = "black" )
segments(x0 = 2.7, y0 =  M_larch, x1 = 3.3, y1= M_larch, lwd = 2, col = "black" )
mtext(c("a", "ab", "b"),side=1,at=c(1,2,3),line = -7, font = 1)


### plot n effectors
#means from transformaton 
M_red<- mean(Red_df_new$n_Effectors)
M_white<- mean(White_df_new$n_Effectors)
M_larch<- mean(Larch_df_new$n_Effectors)

par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))
stripchart(all_df$n_Effectors ~ all_df$group,
           vertical = TRUE,
           method = "jitter", jitter = 0.2, 
           pch = 16, 
           col = c("#5676A166","#B0913666", "#B3675766"),
           bg = rep(c("#5676A166","#B0913666", "#B3675766")),
           cex.axis = 0.7,
           ylim=c(0,250), 
           ylab = expression(paste("Effectors")), 
           axes = FALSE, 
           cex = 1.3)
box()
axis(2)
mtext(text = c("Red", "White", "Larch"),side=1,at=c(1,2,3),line = 1, font = 1)
segments(x0 = .7, y0 =  M_red, x1 = 1.3, y1= M_red, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  M_white, x1 = 2.3, y1= M_white, lwd = 2, col = "black" )
segments(x0 = 2.7, y0 =  M_larch, x1 = 3.3, y1= M_larch, lwd = 2, col = "black" )
mtext(c("a", "a", "a"),side=1,at=c(1,2,3),line = -7, font = 1)


### plot %SSPs out of total proteins
#means from transformaton 
M_red<- mean(Red_df_new$percent_SSPs_out_of_total_genes)
M_white<- mean(White_df_new$percent_SSPs_out_of_total_genes)
M_larch<- mean(Larch_df_new$percent_SSPs_out_of_total_genes)

par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))
stripchart(all_df$percent_SSPs_out_of_total_genes ~ all_df$group,
           vertical = TRUE,
           method = "jitter", jitter = 0.2, 
           pch = 16, 
           col = c("#5676A166","#B0913666", "#B3675766"),
           bg = rep(c("#5676A166","#B0913666", "#B3675766")),
           cex.axis = 0.7,
           ylim=c(0,5), 
           ylab = expression(paste("% SSPs out of all proteins")), 
           axes = FALSE, 
           cex = 1.3)
box()
axis(2)
mtext(text = c("Red", "White", "Larch"),side=1,at=c(1,2,3),line = 1, font = 1)
segments(x0 = .7, y0 =  M_red, x1 = 1.3, y1= M_red, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  M_white, x1 = 2.3, y1= M_white, lwd = 2, col = "black" )
segments(x0 = 2.7, y0 =  M_larch, x1 = 3.3, y1= M_larch, lwd = 2, col = "black" )
mtext(c("a", "a", "a"),side=1,at=c(1,2,3),line = -7, font = 1)


### plot %SSSPs out of SSPs
#means from transformaton 
M_red<- mean(Red_df_new$percent_SSSPs_out_of_SSPs)
M_white<- mean(White_df_new$percent_SSSPs_out_of_SSPs)
M_larch<- mean(Larch_df_new$percent_SSSPs_out_of_SSPs)

par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))
stripchart(all_df$percent_SSSPs_out_of_SSPs ~ all_df$group,
           vertical = TRUE,
           method = "jitter", jitter = 0.2, 
           pch = 16, 
           col = c("#5676A166","#B0913666", "#B3675766"),
           bg = rep(c("#5676A166","#B0913666", "#B3675766")),
           cex.axis = 0.7,
           ylim=c(0,50), 
           ylab = expression(paste("% SSSPs out of SSPs")), 
           axes = FALSE, 
           cex = 1.3)
box()
axis(2)
mtext(text = c("Red", "White", "Larch"),side=1,at=c(1,2,3),line = 1, font = 1)
segments(x0 = .7, y0 =  M_red, x1 = 1.3, y1= M_red, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  M_white, x1 = 2.3, y1= M_white, lwd = 2, col = "black" )
segments(x0 = 2.7, y0 =  M_larch, x1 = 3.3, y1= M_larch, lwd = 2, col = "black" )
mtext(c("a", "ab", "b"),side=1,at=c(1,2,3),line = -7, font = 1)


### plot % Effectors out of SSPs
#means from transformaton 
M_red<- mean(Red_df_new$percent_effectors_out_of_SSPs_inv_square)
M_white<- mean(White_df_new$percent_effectors_out_of_SSPs_inv_square)
M_larch<- mean(Larch_df_new$percent_effectors_out_of_SSPs_inv_square)

par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))
stripchart(all_df$percent_effectors_out_of_SSPs_inv_square ~ all_df$group,
           vertical = TRUE,
           method = "jitter", jitter = 0.2, 
           pch = 16, 
           col = c("#5676A166","#B0913666", "#B3675766"),
           bg = rep(c("#5676A166","#B0913666", "#B3675766")),
           cex.axis = 0.7,
           ylim=c(0,.0018), 
           ylab = expression(paste("% Effectors out of SSPs"^-2)), 
           axes = FALSE, 
           cex = 1.3)
box()
axis(2)
mtext(text = c("Red", "White", "Larch"),side=1,at=c(1,2,3),line = 1, font = 1)
segments(x0 = .7, y0 =  M_red, x1 = 1.3, y1= M_red, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  M_white, x1 = 2.3, y1= M_white, lwd = 2, col = "black" )
segments(x0 = 2.7, y0 =  M_larch, x1 = 3.3, y1= M_larch, lwd = 2, col = "black" )
mtext(c("a", "a", "a"),side=1,at=c(1,2,3),line = -7, font = 1)
