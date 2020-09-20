#load libraries
library("seqinr")
library("data.table")
library("stringr")
library("dplyr")
library("car")
library("FSA")
library("rcompanion")
library("coin")
library("multcompView")
library("multcomp")

#to load coin
if(!require(coin)){install.packages("coin")}
if(!require(FSA)){install.packages("FSA")}
if(!require(rcompanion)){install.packages("rcompanion")}
if(!require(multcompView)){install.packages("multcompView")}

#set options
options(stringsAsFactors = FALSE)

#set wd
setwd("~/Desktop/Manuscript_NP_resubmission/Code/")

#read in data 
input<-read.delim("SSP_prot_genomesize_totals.csv", header = TRUE, sep = ",", fill = TRUE, strip.white = TRUE)

###S vs O

#separate Suillus from Other
input_Suillus<- data.frame(input[input$treatment == "S", ])
input_Other<- data.frame(input[input$treatment == "O", ])

##Stats for genome size 
Suillus_genome_size_mean<- mean(input_Suillus$genome_size)
Suillus_genome_size_mean
Other_genome_size_mean<- mean(input_Other$genome_size)
Other_genome_size_mean

#standard error funciton
std <- function(x) sd(x)/sqrt(length(x))
Suillus_genome_size_SE<- std(input_Suillus$genome_size)
Suillus_genome_size_SE
Other_genome_size_SE<- std(input_Other$genome_size)
Other_genome_size_SE

##is it significant?

#test normality  - if shapiro p-val is < 0.05, its not normal 
shapiro.test(input_Suillus$genome_size)
#not normal
shapiro.test(input_Other$genome_size)
#not normal

#test for equal variance  - if p-val is < 0.05 var is sig different
var.test(input_Suillus$genome_size, input_Other$genome_size)
#variance is not sig. dif, use = var in t-test

#transform to improve normal
input_Suillus$log_genome_size<- log(input_Suillus$genome_size)
input_Other$log_genome_size<- log(input_Other$genome_size)
#re-test normality 
shapiro.test(input_Suillus$log_genome_size)
#normal now!
shapiro.test(input_Other$log_genome_size)
#not quite, but much imporved

#re-test for equal variance  - if p-val is < 0.05 var is sig different
var.test(input_Suillus$log_genome_size, input_Other$log_genome_size)
#variance is not sig. dif, use = var in t-test

#t-test 
t.test(input_Suillus$log_genome_size, input_Other$log_genome_size, var.equal = TRUE)
#not significantly different

##Stats for protein number
Suillus_n_proteins_mean<- mean(input_Suillus$n_proteins)
Suillus_n_proteins_mean
Other_n_proteins_mean<- mean(input_Other$n_proteins)
Other_n_proteins_mean

#standard error function
Suillus_n_proteins_SE<- std(input_Suillus$n_proteins)
Suillus_n_proteins_SE
Other_n_proteins_SE<- std(input_Other$n_proteins)
Other_n_proteins_SE

##is it significant?

#test normality  - if shapiro p-val is < 0.05, its not normal 
shapiro.test(input_Suillus$n_proteins)
#not normal
shapiro.test(input_Other$n_proteins)
#normal

#test for equal variance  - if p-val is < 0.05 var is sig different
var.test(input_Suillus$n_proteins, input_Other$n_proteins)
#variance is sig. dif, don't use = var in t-test

#transform to improve normal
input_Suillus$log_n_proteins<- log(input_Suillus$n_proteins)
input_Other$log_n_proteins<- log(input_Other$n_proteins)
#re-test normality 
shapiro.test(input_Suillus$log_n_proteins)
#normal now!
shapiro.test(input_Other$log_n_proteins)
#still normal

#re-test for equal variance  - if p-val is < 0.05 var is sig different
var.test(input_Suillus$log_n_proteins, input_Other$log_n_proteins)
#variance is sig. dif, don't use = var in t-test

#t-test 
t.test(input_Suillus$log_n_proteins, input_Other$log_n_proteins, var.equal = FALSE)
#not significantly different


##Stats for n SSPs
Suillus_n_SSPs_mean<- mean(input_Suillus$n_SSPs)
Suillus_n_SSPs_mean
Other_n_SSPs_mean<- mean(input_Other$n_SSPs)
Other_n_SSPs_mean

#standard error function
Suillus_n_SSPs_SE<- std(input_Suillus$n_SSPs)
Suillus_n_SSPs_SE
Other_n_SSPs_SE<- std(input_Other$n_SSPs)
Other_n_SSPs_SE

##is it significant?

#test normality  - if shapiro p-val is < 0.05, its not normal 
shapiro.test(input_Suillus$n_SSPs)
#not normal
shapiro.test(input_Other$n_SSPs)
#normal

#test for equal variance  - if p-val is < 0.05 var is sig different
var.test(input_Suillus$n_SSPs, input_Other$n_SSPs)
#variance is sig. dif, don't use = var in t-test

#transform to improve normal
input_Suillus$log_n_SSPs<- log(input_Suillus$n_SSPs)
input_Other$log_n_SSPs<- log(input_Other$n_SSPs)
#re-test normality 
shapiro.test(input_Suillus$log_n_SSPs)
#normal now!
shapiro.test(input_Other$log_n_SSPs)
#still normal

#re-test for equal variance  - if p-val is < 0.05 var is sig different
var.test(input_Suillus$log_n_SSPs, input_Other$log_n_SSPs)
#variance is sig. dif, don't use = var in t-test

#t-test 
t.test(input_Suillus$log_n_SSPs, input_Other$log_n_SSPs, var.equal = FALSE)
#not significantly different


##Stats for n SSSPs
Suillus_n_SSSPs_mean<- mean(input_Suillus$n_SSSPs)
Suillus_n_SSSPs_mean
Other_n_SSSPs_mean<- mean(input_Other$n_SSSPs)
Other_n_SSSPs_mean

#standard error function
Suillus_n_SSSPs_SE<- std(input_Suillus$n_SSSPs)
Suillus_n_SSSPs_SE
Other_n_SSSPs_SE<- std(input_Other$n_SSSPs)
Other_n_SSSPs_SE

##is it significant?

#test normality  - if shapiro p-val is < 0.05, its not normal 
shapiro.test(input_Suillus$n_SSSPs)
#normal
shapiro.test(input_Other$n_SSSPs)
#normal

#test for equal variance  - if p-val is < 0.05 var is sig different
var.test(input_Suillus$n_SSSPs, input_Other$n_SSSPs)
#variance is sig. dif, don't use = var in t-test

#t-test 
t.test(input_Suillus$n_SSSPs, input_Other$n_SSSPs, var.equal = FALSE)
#not significantly different


##Stats for percent SSPs out of proteins
Suillus_percent_SSP_out_of_proteins_mean<- mean(input_Suillus$percent_SSP_out_of_proteins)
Suillus_percent_SSP_out_of_proteins_mean
Other_percent_SSP_out_of_proteins_mean<- mean(input_Other$percent_SSP_out_of_proteins)
Other_percent_SSP_out_of_proteins_mean

#standard error function
Suillus_percent_SSP_out_of_proteins_SE<- std(input_Suillus$percent_SSP_out_of_proteins)
Suillus_percent_SSP_out_of_proteins_SE
Other_percent_SSP_out_of_proteins_SE<- std(input_Other$percent_SSP_out_of_proteins)
Other_percent_SSP_out_of_proteins_SE

##is it significant?

#test normality  - if shapiro p-val is < 0.05, its not normal 
shapiro.test(input_Suillus$percent_SSP_out_of_proteins)
#normal
shapiro.test(input_Other$percent_SSP_out_of_proteins)
#normal

#test for equal variance  - if p-val is < 0.05 var is sig different
var.test(input_Suillus$percent_SSP_out_of_proteins, input_Other$percent_SSP_out_of_proteins)
#variance is not sig. dif, use = var in t-test

#t-test 
t.test(input_Suillus$percent_SSP_out_of_proteins, input_Other$percent_SSP_out_of_proteins, var.equal = TRUE)
#not significantly different

##Stats for percent SSSPs out of SSPs
Suillus_percent_SSSPs_out_of_SSPs_mean<- mean(input_Suillus$percent_SSSPs_out_of_SSPs)
Suillus_percent_SSSPs_out_of_SSPs_mean
Other_percent_SSSPs_out_of_SSPs_mean<- mean(input_Other$percent_SSSPs_out_of_SSPs)
Other_percent_SSSPs_out_of_SSPs_mean

#standard error function
Suillus_percent_SSSPs_out_of_SSPs_SE<- std(input_Suillus$percent_SSSPs_out_of_SSPs)
Suillus_percent_SSSPs_out_of_SSPs_SE
Other_percent_SSSPs_out_of_SSPs_SE<- std(input_Other$percent_SSSPs_out_of_SSPs)
Other_percent_SSSPs_out_of_SSPs_SE

##is it significant?

#test normality  - if shapiro p-val is < 0.05, its not normal 
shapiro.test(input_Suillus$percent_SSSPs_out_of_SSPs)
#normal
shapiro.test(input_Other$percent_SSSPs_out_of_SSPs)
#normal

#test for equal variance  - if p-val is < 0.05 var is sig different
var.test(input_Suillus$percent_SSSPs_out_of_SSPs, input_Other$percent_SSSPs_out_of_SSPs)
#variance is sig. dif, don't use = var in t-test

#t-test 
t.test(input_Suillus$percent_SSSPs_out_of_SSPs, input_Other$percent_SSSPs_out_of_SSPs, var.equal = FALSE)
#not significantly different


###make figures for S vs O 
#order
input$treatment <- factor(input$treatment, levels = c("S", "O"))
#set par
par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0), las=1)

#n SSPs
stripchart(input$n_SSPs ~ input$treatment,
           vertical = TRUE,
           method = "jitter", jitter = 0.2, 
           pch = 16, 
           col = c("#40595266","#9B905B66"),
           bg = rep(c("#40595266","#9B905B66")),
           cex.axis = 0.7,
           ylim=c(0,700), 
           ylab = "n_SSPs", 
           axes = FALSE, 
           cex = 1.3)
box()
axis(2)
mtext(text = c("Suillus", "Other ECM"),side=1,at=c(1,2),line = 1, font = 3)
segments(x0 = .7, y0 =  Suillus_n_SSPs_mean, x1 = 1.3, y1= Suillus_n_SSPs_mean, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  Other_n_SSPs_mean, x1 = 2.3, y1= Other_n_SSPs_mean, lwd = 2, col = "black" )
mtext(c("a", "a"),side=1,at=c(1,2),line = -9, font = 3)

#n_SSSPs
stripchart(input$n_SSSPs ~ input$treatment,
           vertical = TRUE,
           method = "jitter", jitter = 0.2, 
           pch = 16, 
           col = c("#40595266","#9B905B66"),
           bg = rep(c("#40595266","#9B905B66")),
           cex.axis = 0.7,
           ylim=c(0,600), 
           ylab = "n_SSSPs", 
           axes = FALSE, 
           cex = 1.3)
box()
axis(2)
mtext(text = c("Suillus", "Other ECM"),side=1,at=c(1,2),line = 1, font = 3)
segments(x0 = .7, y0 =  Suillus_n_SSSPs_mean, x1 = 1.3, y1= Suillus_n_SSSPs_mean, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  Other_n_SSSPs_mean, x1 = 2.3, y1= Other_n_SSSPs_mean, lwd = 2, col = "black" )
mtext(c("a", "b"),side=1,at=c(1,2),line = -9, font = 3)


#percent_SSSPs_out_of_SSPs
stripchart(input$percent_SSSPs_out_of_SSPs ~ input$treatment,
           vertical = TRUE,
           method = "jitter", jitter = 0.2, 
           pch = 16, 
           col = c("#40595266","#9B905B66"),
           bg = rep(c("#40595266","#9B905B66")),
           cex.axis = 0.7,
           ylim=c(0,90), 
           ylab = "percent_SSSPs_out_of_SSPs", 
           axes = FALSE, 
           cex = 1.3)
box()
axis(2)
mtext(text = c("Suillus", "Other ECM"),side=1,at=c(1,2),line = 1, font = 3)
segments(x0 = .7, y0 =  Suillus_percent_SSSPs_out_of_SSPs_mean, x1 = 1.3, y1= Suillus_percent_SSSPs_out_of_SSPs_mean, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  Other_percent_SSSPs_out_of_SSPs_mean, x1 = 2.3, y1= Other_percent_SSSPs_out_of_SSPs_mean, lwd = 2, col = "black" )
mtext(c("a", "b"),side=1,at=c(1,2),line = -9, font = 3)




###inter-genus comparison###

#split the two data categories for the three groups and get stats for each 
Red_df<- input[input$host=="R",]
White_df<- input[input$host=="W",]
Larch_df<- input[input$host=="L",]

#check that the separation worked
nrow(Red_df)
nrow(White_df)
nrow(Larch_df)
#yep

###genome size 
#means
Red_genome_size_mean<- mean(Red_df$genome_size)
White_genome_size_mean<- mean(White_df$genome_size)
Larch_genome_size_mean<- mean(Larch_df$genome_size)

#SE
Red_genome_size_SE<- std(Red_df$genome_size) 
Red_genome_size_SE
White_genome_size_SE<- std(White_df$genome_size) 
White_genome_size_SE
Larch_genome_size_SE<- std(Larch_df$genome_size)  
Larch_genome_size_SE

#bind all 
all_df<- rbind(Red_df, White_df, Larch_df)
View(all_df)
#build model
genome_size_model <- lm(all_df$genome_size ~ all_df$host)

#test normality and var
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(genome_size_model) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
boxCox(genome_size_model, lambda = seq(-5, 1, 1/10)) #what does boxcox suggest? log transformation works. 

#log transform to improve normality?
all_df$genome_size_log<- log(all_df$genome_size)

#see if the transformation helped
genome_size_model2 <- lm(all_df$genome_size_log ~ all_df$host)
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(genome_size_model2) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
#not a lot, try -2 ==(1/x^2) - inverse square transformation 

all_df$genome_size_inv_suqare<- all_df$genome_size^(-2)

genome_size_model3 <- lm(all_df$genome_size_inv_suqare ~ all_df$host)
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
Red_prot_SE<- std(Red_df$n_proteins)  
Red_prot_SE
White_prot_SE<- std(White_df$n_proteins) 
White_prot_SE
Larch_prot_SE<- std(Larch_df$n_proteins)  
Larch_prot_SE

#build model
prot_model <- lm(all_df$n_proteins ~ all_df$host)

#vis normality and variance
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(prot_model) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
boxCox(prot_model) #what does boxcox suggest? log transformation works. 
powerTransform(prot_model) #how about powerTransform? #log should be fine for this.

#log transform to improve normality?
all_df$n_proteins_log<- log(all_df$n_proteins)

#see if the transformation helped
prot_model2 <- lm(all_df$n_proteins_log ~ all_df$host)
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(prot_model2) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
#well, that didn't help...

#try the recomended transformation
prot_model3 <- lm(all_df$n_proteins^-.5 ~ all_df$host)
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(prot_model3) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
#worse! go back. do not pass go. do not colelct $100. Use original data. 

#run anova
summary(aov(prot_model))
#not significant

#make df
r_GS<- data.frame(Red_df$n_proteins, rep("Red"))
colnames(r_GS)<- c("GS", "host")
w_GS<- data.frame(White_df$n_proteins, rep("White"))
colnames(w_GS)<- c("GS", "host")
l_GS<- data.frame(Larch_df$n_proteins, rep("Larch"))
colnames(l_GS)<- c("GS", "host")
#bind them
GS_df<- rbind(r_GS, w_GS, l_GS)
#set factors
#set order of boxes
GS_df$host <- factor(GS_df$host,levels = c('Red','White', 'Larch'),ordered = TRUE)

#run independence test (this is like an ANOVA)
independence_test(GS ~ host, 
                  data = GS_df)


####SSPs
#means
Red_ssp_mean<- mean(Red_df$n_SSPs)
White_ssp_mean<- mean(White_df$n_SSPs)
Larch_ssp_mean<- mean(Larch_df$n_SSPs)

#SD and SE
Red_ssp_SE<- std(Red_df$n_SSPs) 
Red_ssp_SE
White_ssp_SE<- std(White_df$n_SSPs)  
White_ssp_SE
Larch_ssp_SE<- std(Larch_df$n_SSPs)  
Larch_ssp_SE

#build model
ssp_model <- lm(all_df$n_SSPs ~ all_df$host)

#vis normality and variance
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(ssp_model) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
boxCox(ssp_model, lambda = seq(-6,2,2)) #what does boxcox suggest?  
powerTransform(ssp_model) #how about powerTransform? 
# log to ^-2 

all_df$n_SSPs_log<- log(all_df$n_SSPs)

ssp_model2 <- lm(all_df$n_SSPs_log ~ all_df$host)

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
colnames(r_GS)<- c("GS", "host")
w_GS<- data.frame(White_df$n_SSPs, rep("White"))
colnames(w_GS)<- c("GS", "host")
l_GS<- data.frame(Larch_df$n_SSPs, rep("Larch"))
colnames(l_GS)<- c("GS", "host")
#bind them
GS_df<- rbind(r_GS, w_GS, l_GS)
#set factors
#set order of boxes
GS_df$host <- factor(GS_df$host,levels = c('Red','White', 'Larch'),ordered = TRUE)

#run independence test (this is like an ANOVA)
independence_test(GS ~ host, 
                  data = GS_df)


####SSSPs
#means
Red_sssp_mean<- mean(Red_df$n_SSSPs)
White_sssp_mean<- mean(White_df$n_SSSPs)
Larch_sssp_mean<- mean(Larch_df$n_SSSPs)

#SD and SE
Red_sssp_SE<- std(Red_df$n_SSSPs)  
Red_sssp_SE
White_sssp_SE<- std(White_df$n_SSSPs)
White_sssp_SE
Larch_sssp_SE<- std(Larch_df$n_SSSPs) 
Larch_sssp_SE

#build model
sssp_model <- lm(all_df$n_SSSPs ~ all_df$host)

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
r_GS<- data.frame(Red_df$n_SSSPs, rep("Red"))
colnames(r_GS)<- c("GS", "host")
w_GS<- data.frame(White_df$n_SSSPs, rep("White"))
colnames(w_GS)<- c("GS", "host")
l_GS<- data.frame(Larch_df$n_SSSPs, rep("Larch"))
colnames(l_GS)<- c("GS", "host")
#bind them
GS_df<- rbind(r_GS, w_GS, l_GS)
#set factors
#set order of boxes
GS_df$host <- factor(GS_df$host,levels = c('Red','White', 'Larch'),ordered = TRUE)

#run independence test (this is like an ANOVA)
independence_test(GS ~ host, 
                  data = GS_df)
#significant!

#run tukey equivalent to see whats significant
#take a look at the medians 
boxplot(GS ~ host,
        data = GS_df)

#order them by median (highest to lowest)
GS_df$host = factor(GS_df$host , 
                    levels = c("Larch", "White", "Red"))


headtail(GS_df)

### Pairwise tests
PT = pairwisePermutationTest(GS ~ host,
                             data = GS_df,
                             method="fdr")
PT


#%SSPs out of all prot.
#means
Red_SSP_per_all_prot_mean<- mean(Red_df$percent_SSP_out_of_proteins)
White_SSP_per_all_prot_mean<- mean(White_df$percent_SSP_out_of_proteins)
Larch_SSP_per_all_prot_mean<- mean(Larch_df$percent_SSP_out_of_proteins)

#SD and SE
Red_SSP_per_all_prot_SE<- std(Red_df$percent_SSP_out_of_proteins) 
Red_SSP_per_all_prot_SE
White_SSP_per_all_prot_SE<- std(White_df$percent_SSP_out_of_proteins)  
White_SSP_per_all_prot_SE
Larch_SSP_per_all_prot_SE<- std(Larch_df$percent_SSP_out_of_proteins) 
Larch_SSP_per_all_prot_SE

#build model
SSP_per_all_prot_model <- lm(all_df$percent_SSP_out_of_proteins ~ all_df$host)

#vis normality and variance
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(SSP_per_all_prot_model) #data is non normal, variance looks ok-ish
par(mfrow = c(1,1)) #back to one image per page
boxCox(SSP_per_all_prot_model, lambda = seq(-6,2,2)) #what does boxcox suggest?  
powerTransform(SSP_per_all_prot_model) #how about powerTransform? # ^-2 is what's recommended, but zero is in there, try log.

#transform
all_df$percent_SSP_out_of_proteins_log<- log(all_df$percent_SSP_out_of_proteins)

#make new model
SSP_per_all_prot_model2 <- lm(all_df$percent_SSP_out_of_proteins_log ~ all_df$host)

#see if the transformation helped
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(SSP_per_all_prot_model2) #data looks better
par(mfrow = c(1,1)) #back to one image per page
#looks about the same, just use the original


#run anova
summary(aov(SSP_per_all_prot_model))
#not significant

#randomization test make df
r_GS<- data.frame(Red_df$percent_SSP_out_of_proteins, rep("Red"))
colnames(r_GS)<- c("GS", "host")
w_GS<- data.frame(White_df$percent_SSP_out_of_proteins, rep("White"))
colnames(w_GS)<- c("GS", "host")
l_GS<- data.frame(Larch_df$percent_SSP_out_of_proteins, rep("Larch"))
colnames(l_GS)<- c("GS", "host")
#bind them
GS_df<- rbind(r_GS, w_GS, l_GS)
#set factors
#set order of boxes
GS_df$host <- factor(GS_df$host,levels = c('Red','White', 'Larch'),ordered = TRUE)

#run independence test (this is like an ANOVA)
independence_test(GS ~ host, 
                  data = GS_df)
#marganaly significant 

#whats different?
#run tukey equivelent to see whats significant
#take a look at the medians 
boxplot(GS ~ host,
        data = GS_df)

#order them by median (highest to lowest)
GS_df$host = factor(GS_df$host , 
                    levels = c("Larch", "White", "Red"))


headtail(GS_df)

### Pairwise tests
PT = pairwisePermutationTest(GS ~ host,
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
Red_SSSP_per_SSP_SE<- std(Red_df$percent_SSSPs_out_of_SSPs) 
Red_SSSP_per_SSP_SE
White_SSSP_per_SSP_SE<- std(White_df$percent_SSSPs_out_of_SSPs) 
White_SSSP_per_SSP_SE
Larch_SSSP_per_SSP_SE<- std(Larch_df$percent_SSSPs_out_of_SSPs) 
Larch_SSSP_per_SSP_SE

#build model
SSSP_per_SSP_model <- lm(all_df$percent_SSSPs_out_of_SSPs ~ all_df$host)

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
colnames(r_GS)<- c("GS", "host")
w_GS<- data.frame(White_df$percent_SSSPs_out_of_SSPs, rep("White"))
colnames(w_GS)<- c("GS", "host")
l_GS<- data.frame(Larch_df$percent_SSSPs_out_of_SSPs, rep("Larch"))
colnames(l_GS)<- c("GS", "host")
#bind them
GS_df<- rbind(r_GS, w_GS, l_GS)
#set factors
#set order of boxes
GS_df$host <- factor(GS_df$host,levels = c('Red','White', 'Larch'),ordered = TRUE)

#run independence test (this is like an ANOVA)
independence_test(GS ~ host, 
                  data = GS_df)

#significant!


#run tukey equivelent to see whats significant
#take a look at the medians 
boxplot(GS ~ host,
        data = GS_df)

#order them by median (highest to lowest)
GS_df$host = factor(GS_df$host , 
                    levels = c("Larch", "White", "Red"))


headtail(GS_df)

### Pairwise tests
PT = pairwisePermutationTest(GS ~ host,
                             data = GS_df,
                             method="fdr")
PT


######## plots start here ########

#order red -> white -> larch 
all_df$host <- factor(all_df$host, levels = c("R", "W", "L"))

#set par
par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))

#plot of n_SSPs
stripchart(all_df$n_SSPs ~ all_df$host,
           vertical = TRUE,
           method = "jitter", jitter = 0.2, 
           pch = 16, 
           col = c("#5676A166","#B0913666", "#B3675766"),
           bg = rep(c("#5676A166","#B0913666", "#B3675766")),
           cex.axis = 0.7,
           ylim=c(0,600), 
           ylab = expression(paste("SSPs")), 
           axes = FALSE, 
           cex = 1.3)
box()
axis(2)
mtext(text = c("Red", "White", "Larch"),side=1,at=c(1,2,3),line = 1, font = 1)
segments(x0 = .7, y0 =  Red_ssp_mean, x1 = 1.3, y1= Red_ssp_mean, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  White_ssp_mean, x1 = 2.3, y1= White_ssp_mean, lwd = 2, col = "black" )
segments(x0 = 2.7, y0 =  Larch_ssp_mean, x1 = 3.3, y1= Larch_ssp_mean, lwd = 2, col = "black" )
mtext(c("a", "a", "a"),side=1,at=c(1,2,3),line = -9, font = 1)


###plot of n SSSPs

#means from transformaton 
stripchart(all_df$n_SSSPs ~ all_df$host,
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
segments(x0 = .7, y0 =  Red_sssp_mean, x1 = 1.3, y1= Red_sssp_mean, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  White_sssp_mean, x1 = 2.3, y1= White_sssp_mean, lwd = 2, col = "black" )
segments(x0 = 2.7, y0 =  Larch_sssp_mean, x1 = 3.3, y1= Larch_sssp_mean, lwd = 2, col = "black" )
mtext(c("a", "ab", "b"),side=1,at=c(1,2,3),line = -9, font = 1)


### plot %SSSPs out of SSPs
#means from transformaton 
stripchart(all_df$percent_SSSPs_out_of_SSPs ~ all_df$host,
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
segments(x0 = .7, y0 =  Red_SSSP_per_SSP_mean, x1 = 1.3, y1= Red_SSSP_per_SSP_mean, lwd = 2, col = "black" )
segments(x0 = 1.7, y0 =  White_SSSP_per_SSP_mean, x1 = 2.3, y1= White_SSSP_per_SSP_mean, lwd = 2, col = "black" )
segments(x0 = 2.7, y0 =  Larch_SSSP_per_SSP_mean, x1 = 3.3, y1= Larch_SSSP_per_SSP_mean, lwd = 2, col = "black" )
mtext(c("a", "ab", "b"),side=1,at=c(1,2,3),line = -9, font = 1)

