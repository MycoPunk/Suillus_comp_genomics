#load libraries
library("seqinr")
library("data.table")
library("stringr")
library("car")
library("coin")
library("multcompView")
library("multcomp")
library("mvtnorm")

#load coin
#if(!require(coin)){install.packages("coin")}
#if(!require(FSA)){install.packages("FSA")}
#if(!require(rcompanion)){install.packages("rcompanion")}
#if(!require(multcompView)){install.packages("multcompView")}

#install.packages("mvtnorm", repos="http://R-Forge.R-project.org")

options(stringsAsFactors = TRUE)
setwd("~/Desktop/Project_Suillus_comp_genomics/R")

#read it back in 
SSPs_coded_within_Suillus<- read.csv("SSPs_coded_within_Suillus.csv")


#split the two data categories for the three groups and get stats for each 
Red_df<- SSPs_coded_within_Suillus[SSPs_coded_within_Suillus$group=="a_Red",]
White_df<- SSPs_coded_within_Suillus[SSPs_coded_within_Suillus$group=="b_White",]
Larch_df<- SSPs_coded_within_Suillus[SSPs_coded_within_Suillus$group=="c_Larch",]

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
boxCox(genome_size_model, lambda = seq(-4, 1, 1/10)) #what does boxcox suggest? log transformation works. 

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

####get p val with a randomization test

#make df
r_GS<- data.frame(Red_df$genome_size, rep("Red"))
colnames(r_GS)<- c("GS", "group")
w_GS<- data.frame(White_df$genome_size, rep("White"))
colnames(w_GS)<- c("GS", "group")
l_GS<- data.frame(Larch_df$genome_size, rep("Larch"))
colnames(l_GS)<- c("GS", "group")
#bind them
GS_df<- rbind(r_GS, w_GS, l_GS)

#run permutation test with coin for three factors
independence_test(GS ~ group, 
                  data = GS_df)

####proteins
#means
Red_prot_mean<- mean(Red_df$n_proteins_from_gene_cat)
White_prot_mean<- mean(White_df$n_proteins_from_gene_cat)
Larch_prot_mean<- mean(Larch_df$n_proteins_from_gene_cat)

#SD and SE
Red_prot_SD<- sd(Red_df$n_proteins_from_gene_cat)
Red_prot_SE<- std(Red_df$n_proteins_from_gene_cat)  

White_prot_SD<- sd(White_df$n_proteins_from_gene_cat)
White_prot_SE<- std(White_df$n_proteins_from_gene_cat)  

Larch_prot_SD<- sd(Larch_df$n_proteins_from_gene_cat)
Larch_prot_SE<- std(Larch_df$n_proteins_from_gene_cat)  

#build model
prot_model <- lm(all_df$n_proteins_from_gene_cat ~ all_df$group)

#vis normality and variance
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(prot_model) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
boxCox(prot_model) #what does boxcox suggest? log transformation works. 
powerTransform(prot_model) #how about powerTransform? #log should be find for this.

#log transform to improve normality?
all_df$n_proteins_from_gene_cat_log<- log(all_df$n_proteins_from_gene_cat)

#see if the transformation helped
prot_model2 <- lm(all_df$n_proteins_from_gene_cat_log ~ all_df$group)
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(prot_model2) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
#well, that didn't help...

#try the boc recomended transformation
prot_model3 <- lm(all_df$n_proteins_from_gene_cat^-.5 ~ all_df$group)
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(prot_model3) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
#worse! go back. do not pass go. do not colelct $100. Use original data. 

#run anova
summary(aov(prot_model))
#not significant


##permutation test
#make df
r_prot<- data.frame(Red_df$n_proteins_from_gene_cat, rep("Red"))
colnames(r_prot)<- c("prot", "group")
w_prot<- data.frame(White_df$n_proteins_from_gene_cat, rep("White"))
colnames(w_prot)<- c("prot", "group")
l_prot<- data.frame(Larch_df$n_proteins_from_gene_cat, rep("Larch"))
colnames(l_prot)<- c("prot", "group")
#bind them
prot_df<- rbind(r_prot, w_prot, l_prot)

#run permutation test with coin for three factors
independence_test(prot ~ group, 
                  data = prot_df)

####SSPs
#means
Red_ssp_mean<- mean(Red_df$n_SSPs_signalP_TMHMM_lt_300aa)
White_ssp_mean<- mean(White_df$n_SSPs_signalP_TMHMM_lt_300aa)
Larch_ssp_mean<- mean(Larch_df$n_SSPs_signalP_TMHMM_lt_300aa)

#SD and SE
Red_ssp_SD<- sd(Red_df$n_SSPs_signalP_TMHMM_lt_300aa)
Red_ssp_SE<- std(Red_df$n_SSPs_signalP_TMHMM_lt_300aa)  

White_ssp_SD<- sd(White_df$n_SSPs_signalP_TMHMM_lt_300aa)
White_ssp_SE<- std(White_df$n_SSPs_signalP_TMHMM_lt_300aa)  

Larch_ssp_SD<- sd(Larch_df$n_SSPs_signalP_TMHMM_lt_300aa)
Larch_ssp_SE<- std(Larch_df$n_SSPs_signalP_TMHMM_lt_300aa)  

#build model
ssp_model <- lm(all_df$n_SSPs_signalP_TMHMM_lt_300aa ~ all_df$group)

#vis normality and variance
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(ssp_model) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
boxCox(ssp_model, lambda = seq(-6,2,2)) #what does boxcox suggest?  
powerTransform(ssp_model) #how about powerTransform? # ^-2 is what's recommended 

all_df$n_SSPs_signalP_TMHMM_lt_300aa_inv_square<- all_df$n_SSPs_signalP_TMHMM_lt_300aa^(-2)

ssp_model2 <- lm(all_df$n_SSPs_signalP_TMHMM_lt_300aa_inv_square ~ all_df$group)

#see if the transformation helped
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(ssp_model2) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
#looks kinda better

#run anova
summary(aov(ssp_model2))
#not significant

##permutation test
#make df
r_ssp<- data.frame(Red_df$n_SSPs_signalP_TMHMM_lt_300aa, rep("Red"))
colnames(r_ssp)<- c("ssp", "group")
w_ssp<- data.frame(White_df$n_SSPs_signalP_TMHMM_lt_300aa, rep("White"))
colnames(w_ssp)<- c("ssp", "group")
l_ssp<- data.frame(Larch_df$n_SSPs_signalP_TMHMM_lt_300aa, rep("Larch"))
colnames(l_ssp)<- c("ssp", "group")
#bind them
ssp_df<- rbind(r_ssp, w_ssp, l_ssp)

#run permutation test with coin for three factors
independence_test(ssp ~ group, 
                  data = ssp_df)

####SSSPs
#means
Red_ssp_mean<- mean(Red_df$n_SSSPs)
White_ssp_mean<- mean(White_df$n_SSSPs)
Larch_ssp_mean<- mean(Larch_df$n_SSSPs)

#SD and SE
Red_sssp_SD<- sd(Red_df$n_SSSPs)
Red_sssp_SE<- std(Red_df$n_SSSPs)  

White_sssp_SD<- sd(White_df$n_SSSPs)
White_sssp_SE<- std(White_df$n_SSSPs)  

Larch_sssp_SD<- sd(Larch_df$n_SSSPs)
Larch_sssp_SE<- std(Larch_df$n_SSSPs)  

#build model
sssp_model <- lm(all_df$n_SSSPs ~ all_df$group)

#vis normality and variance
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(sssp_model) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
boxCox(sssp_model, lambda = seq(-6,2,2)) #what does boxcox suggest?  
powerTransform(sssp_model) #how about powerTransform? # ^-1 is what's recommended 

#transform
all_df$n_SSSPs_inv_one<- all_df$n_SSSPs^(-1)

#make new model
sssp_model2 <- lm(all_df$n_SSSPs_inv_one ~ all_df$group)

#see if the transformation helped
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(sssp_model2) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
#looks better

#run anova
summary(aov(sssp_model2))
#significant!

#what's different?
sssp_aov<- aov(sssp_model2)
TukeyHSD(sssp_aov)
#red white, and red and larch are sig. different 

##permutation test
#make df
r_sssp<- data.frame(Red_df$n_SSSPs, rep("Red"))
colnames(r_sssp)<- c("sssp", "group")
w_sssp<- data.frame(White_df$n_SSSPs, rep("White"))
colnames(w_sssp)<- c("sssp", "group")
l_sssp<- data.frame(Larch_df$n_SSSPs, rep("Larch"))
colnames(l_sssp)<- c("sssp", "group")
#bind them
sssp_df<- rbind(r_sssp, w_sssp, l_sssp)

#run permutation test with coin for three factors
independence_test(sssp ~ group, 
                  data = sssp_df)

####Effectors
#means
Red_effector_mean<- mean(Red_df$n_effectors_from_EffectorP)
White_effector_mean<- mean(White_df$n_effectors_from_EffectorP)
Larch_effector_mean<- mean(Larch_df$n_effectors_from_EffectorP)

#SD and SE
Red_effector_SD<- sd(Red_df$n_effectors_from_EffectorP)
Red_effector_SE<- std(Red_df$n_effectors_from_EffectorP)  

White_effector_SD<- sd(White_df$n_effectors_from_EffectorP)
White_effector_SE<- std(White_df$n_effectors_from_EffectorP)  

Larch_effector_SD<- sd(Larch_df$n_effectors_from_EffectorP)
Larch_effector_SE<- std(Larch_df$n_effectors_from_EffectorP)  

#build model
effector_model <- lm(all_df$n_effectors_from_EffectorP ~ all_df$group)

#vis normality and variance
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(effector_model) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
boxCox(effector_model, lambda = seq(-6,2,2)) #what does boxcox suggest?  
powerTransform(effector_model) #how about powerTransform? # ^-2 is what's recommended 

#transform
all_df$n_effectors_from_EffectorP_inv_square<- all_df$n_effectors_from_EffectorP^(-2)

#make new model
effector_model2 <- lm(all_df$n_effectors_from_EffectorP_inv_square ~ all_df$group)

#see if the transformation helped
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(effector_model2) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
#looks better

#run anova
summary(aov(effector_model2))
#not significant

##permutation test
#make df
r_eff<- data.frame(Red_df$n_effectors_from_EffectorP, rep("Red"))
colnames(r_eff)<- c("eff", "group")
w_eff<- data.frame(White_df$n_effectors_from_EffectorP, rep("White"))
colnames(w_eff)<- c("eff", "group")
l_eff<- data.frame(Larch_df$n_effectors_from_EffectorP, rep("Larch"))
colnames(l_eff)<- c("eff", "group")
#bind them
eff_df<- rbind(r_eff, w_eff, l_eff)

#run permutation test with coin for three factors
independence_test(eff ~ group, 
                  data = eff_df)



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
powerTransform(SSP_per_all_prot_model) #how about powerTransform? # ^-2 is what's recommended 

#transform
all_df$percent_SSPs_out_of_total_genes_inv_square<- all_df$percent_SSPs_out_of_total_genes^(-2)

#make new model
SSP_per_all_prot_model2 <- lm(all_df$percent_SSPs_out_of_total_genes_inv_square ~ all_df$group)

#see if the transformation helped
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(SSP_per_all_prot_model2) #data looks better
par(mfrow = c(1,1)) #back to one image per page
#looks better

#run anova
summary(aov(SSP_per_all_prot_model2))
#not significant

##permutation test
#make df
r_eff<- data.frame(Red_df$percent_SSPs_out_of_total_genes, rep("Red"))
colnames(r_eff)<- c("eff", "group")
w_eff<- data.frame(White_df$percent_SSPs_out_of_total_genes, rep("White"))
colnames(w_eff)<- c("eff", "group")
l_eff<- data.frame(Larch_df$percent_SSPs_out_of_total_genes, rep("Larch"))
colnames(l_eff)<- c("eff", "group")
#bind them
eff_df<- rbind(r_eff, w_eff, l_eff)

#run permutation test with coin for three factors
independence_test(eff ~ group, 
                  data = eff_df)


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

##permutation test
#make df
r_eff<- data.frame(Red_df$percent_SSSPs_out_of_SSPs, rep("Red"))
colnames(r_eff)<- c("eff", "group")
w_eff<- data.frame(White_df$percent_SSSPs_out_of_SSPs, rep("White"))
colnames(w_eff)<- c("eff", "group")
l_eff<- data.frame(Larch_df$percent_SSSPs_out_of_SSPs, rep("Larch"))
colnames(l_eff)<- c("eff", "group")
#bind them
eff_df<- rbind(r_eff, w_eff, l_eff)

#run permutation test with coin for three factors
independence_test(eff ~ group, 
                  data = eff_df)


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


##permutation test
#make df
r_eff<- data.frame(Red_df$percent_effectors_out_of_SSPs, rep("Red"))
colnames(r_eff)<- c("eff", "group")
w_eff<- data.frame(White_df$percent_effectors_out_of_SSPs, rep("White"))
colnames(w_eff)<- c("eff", "group")
l_eff<- data.frame(Larch_df$percent_effectors_out_of_SSPs, rep("Larch"))
colnames(l_eff)<- c("eff", "group")
#bind them
eff_df<- rbind(r_eff, w_eff, l_eff)

#run permutation test with coin for three factors
independence_test(eff ~ group, 
                  data = eff_df)



######## plots start here ########
#first remove other host associations from the df
groups<- c("a_Red", "b_White", "c_Larch")
all_df<- droplevels(all_df[, all_df$group %in% groups])



#means from transformaton 
Red_df_new<- all_df[all_df$group=="a_Red",]
White_df_new<- all_df[all_df$group=="b_White",]
Larch_df_new<- all_df[all_df$group=="c_Larch",]

M_red<- mean(Red_df_new$n_SSPs_signalP_TMHMM_lt_300aa_inv_square)
M_white<- mean(White_df_new$n_SSPs_signalP_TMHMM_lt_300aa_inv_square)
M_larch<- mean(Larch_df_new$n_SSPs_signalP_TMHMM_lt_300aa_inv_square)


par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))

stripchart(all_df$n_SSPs_signalP_TMHMM_lt_300aa_inv_square ~ all_df$group,
           vertical = TRUE,
           method = "jitter", jitter = 0.2, 
           pch = 16, 
           col = c("#98ABD8", "#566DB4", "#283666"),
           bg = rep(c("#98ABD8", "#566DB4", "#283666")),
           cex.axis = 0.7,
           ylim=c(0,13E-06), 
           ylab = expression(paste("SSPs"^-2)), 
           axes = FALSE, 
           cex = 1.3)
box()
axis(2)
mtext(text = c("Red", "White", "Larch"),side=1,at=c(1,2,3),line = 1, font = 1)
segments(x0 = .7, y0 =  M_red, x1 = 1.3, y1= M_red, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  M_white, x1 = 2.3, y1= M_white, lwd = 2, col = "black" )
segments(x0 = 2.7, y0 =  M_larch, x1 = 3.3, y1= M_larch, lwd = 2, col = "black" )
mtext(c("a", "a", "a"),side=1,at=c(1,2,3),line = -8, font = 1)


###plot of n SSSPs

#means from transformaton 
M_red<- mean(Red_df_new$n_SSSPs_inv_one)
M_white<- mean(White_df_new$n_SSSPs_inv_one)
M_larch<- mean(Larch_df_new$n_SSSPs_inv_one)

par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))
stripchart(all_df$n_SSSPs_inv_one ~ all_df$group,
           vertical = TRUE,
           method = "jitter", jitter = 0.2, 
           pch = 16, 
           col = c("#98ABD8", "#566DB4", "#283666"),
           bg = rep(c("#98ABD8", "#566DB4", "#283666")),
           cex.axis = 0.7,
           ylim=c(0,0.014), 
           ylab = expression(paste("SSSPs"^-1)), 
           axes = FALSE, 
           cex = 1.3)
box()
axis(2)
mtext(text = c("Red", "White", "Larch"),side=1,at=c(1,2,3),line = 1, font = 1)
segments(x0 = .7, y0 =  M_red, x1 = 1.3, y1= M_red, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  M_white, x1 = 2.3, y1= M_white, lwd = 2, col = "black" )
segments(x0 = 2.7, y0 =  M_larch, x1 = 3.3, y1= M_larch, lwd = 2, col = "black" )
mtext(c("a", "ab", "b"),side=1,at=c(1,2,3),line = -8, font = 1)


### plot n effectors
#means from transformaton 
M_red<- mean(Red_df_new$n_effectors_from_EffectorP_inv_square)
M_white<- mean(White_df_new$n_effectors_from_EffectorP_inv_square)
M_larch<- mean(Larch_df_new$n_effectors_from_EffectorP_inv_square)

par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))
stripchart(all_df$n_effectors_from_EffectorP_inv_square ~ all_df$group,
           vertical = TRUE,
           method = "jitter", jitter = 0.2, 
           pch = 16, 
           col = c("#98ABD8", "#566DB4", "#283666"),
           bg = rep(c("#98ABD8", "#566DB4", "#283666")),
           cex.axis = 0.7,
           ylim=c(0,1.6E-04), 
           ylab = expression(paste("Effectors"^-2)), 
           axes = FALSE, 
           cex = 1.3)
box()
axis(2)
mtext(text = c("Red", "White", "Larch"),side=1,at=c(1,2,3),line = 1, font = 1)
segments(x0 = .7, y0 =  M_red, x1 = 1.3, y1= M_red, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  M_white, x1 = 2.3, y1= M_white, lwd = 2, col = "black" )
segments(x0 = 2.7, y0 =  M_larch, x1 = 3.3, y1= M_larch, lwd = 2, col = "black" )
mtext(c("a", "a", "a"),side=1,at=c(1,2,3),line = -8, font = 1)


### plot %SSPs out of total proteins
#means from transformaton 
M_red<- mean(Red_df_new$percent_SSPs_out_of_total_genes_inv_square)
M_white<- mean(White_df_new$percent_SSPs_out_of_total_genes_inv_square)
M_larch<- mean(Larch_df_new$percent_SSPs_out_of_total_genes_inv_square)

par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))
stripchart(all_df$percent_SSPs_out_of_total_genes_inv_square ~ all_df$group,
           vertical = TRUE,
           method = "jitter", jitter = 0.2, 
           pch = 16, 
           col = c("#98ABD8", "#566DB4", "#283666"),
           bg = rep(c("#98ABD8", "#566DB4", "#283666")),
           cex.axis = 0.7,
           ylim=c(0,.4), 
           ylab = expression(paste("% SSPs out of all proteins"^-2)), 
           axes = FALSE, 
           cex = 1.3)
box()
axis(2)
mtext(text = c("Red", "White", "Larch"),side=1,at=c(1,2,3),line = 1, font = 1)
segments(x0 = .7, y0 =  M_red, x1 = 1.3, y1= M_red, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  M_white, x1 = 2.3, y1= M_white, lwd = 2, col = "black" )
segments(x0 = 2.7, y0 =  M_larch, x1 = 3.3, y1= M_larch, lwd = 2, col = "black" )
mtext(c("a", "a", "a"),side=1,at=c(1,2,3),line = -8, font = 1)


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
           col = c("#98ABD8", "#566DB4", "#283666"),
           bg = rep(c("#98ABD8", "#566DB4", "#283666")),
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
mtext(c("a", "ab", "b"),side=1,at=c(1,2,3),line = -8, font = 1)


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
           col = c("#98ABD8", "#566DB4", "#283666"),
           bg = rep(c("#98ABD8", "#566DB4", "#283666")),
           cex.axis = 0.7,
           ylim=c(0,.0016), 
           ylab = expression(paste("% Effectors out of SSPs"^-2)), 
           axes = FALSE, 
           cex = 1.3)
box()
axis(2)
mtext(text = c("Red", "White", "Larch"),side=1,at=c(1,2,3),line = 1, font = 1)
segments(x0 = .7, y0 =  M_red, x1 = 1.3, y1= M_red, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  M_white, x1 = 2.3, y1= M_white, lwd = 2, col = "black" )
segments(x0 = 2.7, y0 =  M_larch, x1 = 3.3, y1= M_larch, lwd = 2, col = "black" )
mtext(c("a", "a", "a"),side=1,at=c(1,2,3),line = -8, font = 1)
