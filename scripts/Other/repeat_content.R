setwd("~/Desktop/Project_Suillus_comp_genomics/R")

repeat_df<-read.delim("repeat_content.csv", header = TRUE, sep = ",", fill = TRUE, strip.white = TRUE)

View(repeat_df)


#is there a relationship between repeat content and genome size? 

#first - remove the influence of repeat regions from the genome size totals 
repeat_df$genome_size_wo_repeats<- repeat_df$genome_size..Mbp. - repeat_df$repeat_content_n_Mbp

##take a look
#here it is without accounting for the repeats in genome size
plot(repeat_df$genome_size..Mbp., repeat_df$repeat_content_n_Mbp)
#here's with accounting for that
plot(repeat_df$genome_size_wo_repeats, repeat_df$repeat_content_n_Mbp)

#is it significant? 
dev.off()
#GRAPH Clade1 novel mutations compared at each time point (minus first T pt.)
genome_size_account_for_repeats_plot<- plot(repeat_df$genome_size_wo_repeats, repeat_df$repeat_content_n_Mbp,
                                 col = "red", 
                                 xlab = "genome size (Mbp)", 
                                 ylab = "repeat content (Mbp)",
                                 main = "accounting for repeats", 
                                 pch = 16)

size_lm<- lm(repeat_df$repeat_content_n_Mbp ~ repeat_df$genome_size_wo_repeats)
sum1<- summary(size_lm)
abline(size_lm)
rsq_clade1_new <- bquote(italic(R)^2 == .(format(summary(size_lm)$adj.r.squared, digits=2)))
text(x = 75, y = 20, labels = rsq_clade1_new)
pval_calde1_new<- round(sum1$coefficients[2,4], digits = 4)
pval_calde1_new <- bquote(italic(P) == .(pval_calde1_new))
text(x = 75, y = 25, labels = pval_calde1_new)
#significant even when accounting for the influence of repeats

#get range of repeat content 
max(repeat_df$repeat_content_n_Mbp)
min(repeat_df$repeat_content_n_Mbp)
