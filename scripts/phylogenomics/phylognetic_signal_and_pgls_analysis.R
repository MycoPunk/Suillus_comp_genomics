#this script generates figures and statistics for phylogenetic signal for trait values for the Suillus comparative genomics project

#set directory
setwd("~/Desktop/Manuscript_NP_resubmission/Code/")

#set libraries
library(phylosignal)
library(adephylo)
library(ape)
library(phylobase)
library(phytools)
library(caper)
library(geiger)


##read in trait data
SSP_traits<- read.table(file ="SSP_prot_genomesize_totals.csv", header = TRUE, sep = ",")
SMC_traits<-read.table(file ="AS5_outputs.csv", header = TRUE, sep = ",")
GPCR_traits<- read.table(file ="GPCR_totals.csv", header = TRUE, sep = ",")
GPCR_traits$total<- rowSums(GPCR_traits[5:ncol(GPCR_traits)],)


#make the .phylo file
#read in the tree
nexus<- tree<- "((((((((((((((((((Suivar1:0.010372,Suitom1:0.006211):0.001161,(Suihi1:0.009804,Suifus1:0.009871):0.001006):0.003248,Suiplo1:0.010465):0.00292,Suidis1:0.010283):0.018954,Suibov1:0.042762):0.011931,(((((Suisubl1:0.010257,Suisu1:0.015162):0.013213,Suiame1:0.020213):0.003498,Suisub1:0.0274):0.003311,Suicot1:0.032812):0.00549,(Suipic1:0.020678,Suidec1:0.01842):0.024794):0.003182):0.0023,(Suigr1:0.01368,Suipla1:0.010526):0.016982):0.002008,((Suilu4:0.010432,Suiocc1:0.016253):0.004298,Suibr2:0.011667):0.022327):0.010157,Suilak1:0.033981):0.008243,Suicli1:0.034744):0.017466,(Suipal1:0.037565,Suiamp1:0.034194):0.010046):0.025448,((((Rhivi1:0.013169,Rhives1:0.018102):0.028885,Rhivul1:0.057818):0.006656,Rhisa1:0.051235):0.007496,Rhitru1:0.079972):0.022267):0.113665,(((Pismi1:0.067099,Pisti1:0.055341):0.081864,Sclci1:0.12731):0.085541,(Paxin1:0.06372,Gyrli1:0.062439):0.081453):0.055982):0.104422,Pilcr1:0.205084):0.029249,(((Lacam2:0.030668,Lacbi2:0.031522):0.145046,Hebcy2:0.199426):0.059435,Amamu1:0.281571):0.070896):0.042081,((Theter1:0.080412,Thega1:0.073159):0.28376,(Ruscom1:0.087985,Rusbre1:0.092679):0.251688):0.035305):0.111335,(Hydru2:0.235556,Cananz1:0.352162):0.263901):0.092299,(Gaumor1:0.166361,Hyssto1:0.187679):0.092299);"
tre.w.names<- read.tree(text = tree)


#make phylogenetic distance matrix
tre.w.names_coph<- cophenetic(tre.w.names)

#isolate trait data (SSP etc.)
#SSPs
SSP<- SSP_traits$n_SSPs
names(SSP) = SSP_traits$JGI_project_code
SSSP<- SSP_traits$n_SSSPs
names(SSSP) = SSP_traits$JGI_project_code
percentSSSPoutofSSP<- SSP_traits$percent_SSSPs_out_of_SSPs
names(percentSSSPoutofSSP) = SSP_traits$JGI_project_code

#order the trait values to match the tree order
SSP<- SSP[row.names(tre.w.names_coph)]
SSSP<- SSSP[row.names(tre.w.names_coph)]
percentSSSPoutofSSP<- percentSSSPoutofSSP[row.names(tre.w.names_coph)]

#for SMCs
SMC_all<- SMC_traits$total_BSG_clusters
names(SMC_all) = SMC_traits$JGI_project_code
#order the trait values to match the tree order
SMC_all<- SMC_all[row.names(tre.w.names_coph)]
#terpenes
#isolate trait data
SMC_terp<- SMC_traits$Terpene
names(SMC_terp) = SMC_traits$JGI_project_code
#order the trait values to match the tree order
SMC_terp<- SMC_terp[row.names(tre.w.names_coph)]
#for nrps-like
#isolate trait data
SMC_NRPS_like<- SMC_traits$NRPS.like
names(SMC_NRPS_like) = SMC_traits$JGI_project_code
#order the trait values to match the tree order
SMC_NRPS_like<- SMC_NRPS_like[row.names(tre.w.names_coph)]

#for gpcrs (overall)
#isolate trait data (SMCs)
GPCRs<- GPCR_traits$total
names(GPCRs) = GPCR_traits$JGI_project_code
#order the trait values to match the tree order
GPCRs<- GPCRs[row.names(tre.w.names_coph)]

#set gorups
groups<- SSP_traits$treatment
names(groups) = SSP_traits$JGI_project_code
#order the trait values to match the tree order
groups<- groups[row.names(tre.w.names_coph)]


###phylosgignal graphs for local autocorrelation
#create a dataframe of the three SSP traits 
SSP_input<- data.frame(SSP = SSP, SSSP = SSSP, percentSSSPoutofSSP = percentSSSPoutofSSP)
#create a dataframe of the two SMC traits 
SMC_input<-  data.frame(NRPS_like = SMC_NRPS_like, SMC_terp = SMC_terp)
#create a dataframe of the two GPCRs traits 
GPCR_input<-  data.frame(GPCRs = GPCRs)

#create .p4d object
SSP_input.p4d_named <- phylo4d(tre.w.names, SSP_input)
SMC_input.p4d_named <- phylo4d(tre.w.names, SMC_input)
GPCR_input.p4d_named <- phylo4d(tre.w.names, GPCR_input)

#compute lipaMoran I -SSP
SSP.lipa <- lipaMoran(SSP_input.p4d_named, reps = 10000)
SSP.lipa.p4d <- lipaMoran(SSP_input.p4d_named, as.p4d = TRUE, reps = 10000)

#set colors and cut offs
points.col<-ifelse(SSP.lipa$p.value < 0.01, "#FCD090", "#93959799")
#plot -SSP
barplot.phylo4d(SSP_input.p4d_named, bar.col=points.col, 
                center = TRUE, 
                scale = FALSE, 
                tree.ladderize = TRUE,
                tip.cex = .5,
                trait.cex = .8,
                trait.bg.col = "white",
                grid.vertical = FALSE,
                grid.horizontal = FALSE,
                bar.lwd = 5,
                rotate.tree = 90)

#compute lipaMoran I -SMC
SMC.lipa <- lipaMoran(SMC_input.p4d_named, reps = 10000)
SMC.lipa.p4d <- lipaMoran(SMC_input.p4d_named, as.p4d = TRUE, reps = 10000)

#plot
points.col<-ifelse(SMC.lipa$p.value < 0.01, "#FCD090", "#93959799")
barplot.phylo4d(SMC_input.p4d_named, bar.col=points.col, 
                center = TRUE, 
                scale = FALSE, 
                tree.ladderize = TRUE,
                tip.cex = .5,
                trait.cex = .8,
                trait.bg.col = "white",
                grid.vertical = FALSE,
                grid.horizontal = FALSE,
                bar.lwd = 5)

#compute lipaMoran I -GPCRs
GPCR.lipa <- lipaMoran(GPCR_input.p4d_named, reps = 10000)
GPCR.lipa.p4d <- lipaMoran(GPCR_input.p4d_named, as.p4d = TRUE, reps = 10000)

#plot
points.col<-ifelse(GPCR.lipa$p.value < 0.01, "#FCD090", "#93959799")
barplot.phylo4d(GPCR_input.p4d_named, bar.col=points.col, 
                center = TRUE, 
                scale = FALSE, 
                tree.ladderize = TRUE,
                tip.cex = .5,
                trait.cex = .8,
                trait.bg.col = "white",
                grid.vertical = FALSE,
                grid.horizontal = FALSE,
                bar.lwd = 3.5)


###compute phylogenetic signal at each internal node using phyloSignalINT()
#ssp
ssplike_internal_node<- phyloSignalINT(SSP_input.p4d_named, trait = names(tipData(SSP_input.p4d_named))[1], method = "K",
                                       reps = 1000, W = NULL)

ssplike_node<- tdata(ssplike_internal_node, "internal")[5]
#make is a named list
ssplike_node_list<- ssplike_node$pvalue.K.SSP
names(ssplike_node_list) <- rownames(ssplike_node)
#reduce significant figures
ssplike_node_list<-round(ssplike_node_list, digits = 4)
#plot
plot(tre.w.names, cex = 0.5)

for (i in 1:length(ssplike_node_list)){
  if (ssplike_node_list[i] < 0.01) 
  {
    nodelabels(ssplike_node_list[[i]],as.integer(names(ssplike_node_list[i])), cex = 0.3, bg="white", frame="c")
  }
}
#conclusion: 52 is significant at 0.005

##SSSPs
sssplike_internal_node<- phyloSignalINT(SSP_input.p4d_named, trait = names(tipData(SSP_input.p4d_named))[2], method = "K",
                                        reps = 1000, W = NULL)

sssplike_node<- tdata(sssplike_internal_node, "internal")[5]
#make is a named list
sssplike_node_list<- sssplike_node$pvalue.K.SSSP
names(sssplike_node_list) <- rownames(sssplike_node)
#reduce significant figures
sssplike_node_list<-round(sssplike_node_list, digits = 4)
#plot
plot(tre.w.names, cex = 0.5)

for (i in 1:length(sssplike_node_list)){
  if (sssplike_node_list[i] < 0.01) 
  {
    nodelabels(sssplike_node_list[[i]],as.integer(names(sssplike_node_list[i])), cex = 0.3, bg="white", frame="c")
  }
}
#conclusion: 53 is significant at 0.007

##percentsssps
percentsssplike_internal_node<- phyloSignalINT(SSP_input.p4d_named, trait = names(tipData(SSP_input.p4d_named))[3], method = "K",
                                               reps = 1000, W = NULL)

percentsssplike_node<- tdata(percentsssplike_internal_node, "internal")[5]
#make is a named list
percentsssplike_node_list<- percentsssplike_node$pvalue.K.percentSSSPoutofSSP
names(percentsssplike_node_list) <- rownames(percentsssplike_node)
#reduce significant figures
percentsssplike_node_list<-round(percentsssplike_node_list, digits = 4)
#plot
plot(tre.w.names, cex = 0.5)

for (i in 1:length(percentsssplike_node_list)){
  if (percentsssplike_node_list[i] < 0.01) 
  {
    nodelabels(percentsssplike_node_list[[i]],as.integer(names(percentsssplike_node_list[i])), cex = 0.3, bg="white", frame="c")
  }
}
#conclusion: 52 is significant at 0.002, 53 is significant at 0.001


##SMCs
#NRPS
nrpslike_internal_node<- phyloSignalINT(SMC_input.p4d_named, trait = names(tipData(SMC_input.p4d_named))[1], method = "K",
                                   reps = 1000, W = NULL)

nrpslike_node<- tdata(nrpslike_internal_node, "internal")[4]
#make is a named list
nrpslike_node_list<- nrpslike_node$pvalue.K.NRPS_like
names(nrpslike_node_list) <- rownames(nrpslike_node)
#reduce significant figures
nrpslike_node_list<-round(nrpslike_node_list, digits = 4)
#plot
plot(tre.w.names, cex = 0.5)

for (i in 1:length(nrpslike_node_list)){
  if (nrpslike_node_list[i] <0.01) 
  {
    nodelabels(nrpslike_node_list[[i]],as.integer(names(nrpslike_node_list[i])), cex = 0.3, bg="white", frame="c")
  }
}
#conclusion: node 52 significant at 0.001

#Terpenes
terpenelike_internal_node<- phyloSignalINT(SMC_input.p4d_named, trait = names(tipData(SMC_input.p4d_named))[2], method = "K",
                                           reps = 1000, W = NULL)

terpenelike_node<- tdata(terpenelike_internal_node, "internal")[4]
#make is a named list
terpenelike_node_list<- terpenelike_node$pvalue.K.SMC_terp
names(terpenelike_node_list) <- rownames(terpenelike_node)
#reduce significant figures
terpenelike_node_list<-round(terpenelike_node_list, digits = 4)
#plot
plot(tre.w.names, cex = 0.5)

for (i in 1:length(terpenelike_node_list)){
  if (terpenelike_node_list[i] <0.01) 
  {
    nodelabels(terpenelike_node_list[[i]],as.integer(names(terpenelike_node_list[i])), cex = 0.3, bg="white", frame="c")
  }
}
#conclusion: node 52 significant at 0.001
            #node 53 significant at 0.001 

#GPCRs
gpcrlike_internal_node<- phyloSignalINT(GPCR_input.p4d_named, trait = names(tipData(GPCR_input.p4d_named))[1], method = "K",
                                        reps = 1000, W = NULL)

gpcrlike_node<- tdata(gpcrlike_internal_node, "internal")[3]
#make is a named list
gpcrlike_node_list<- gpcrlike_node$pvalue.K.GPCRs
names(gpcrlike_node_list) <- rownames(gpcrlike_node)
#reduce significant figures
gpcrlike_node_list<-round(gpcrlike_node_list, digits = 4)
#plot
plot(tre.w.names, cex = 0.5)

for (i in 1:length(gpcrlike_node_list)){
  if (gpcrlike_node_list[i] < 0.01) 
  {
    nodelabels(gpcrlike_node_list[[i]],as.integer(names(gpcrlike_node_list[i])), cex = 0.3, bg="white", frame="c")
  }
}
#conclusion: neither 52 or 53 are significant



###calculate Bloomberg's of the overall phylogeny using phylosignal
#SSP
phyloSignal(SSP_input.p4d_named, methods = "K", reps = 10000)
#SMC
phyloSignal(SMC_input.p4d_named, methods = "K", reps = 10000)
#GPCRs
phyloSignal(GPCR_input.p4d_named, methods = "K", reps = 10000)




#make tree ultimetric for pgls analysis 
library(phangorn) 
force.ultrametric<-function(tree,method=c("nnls","extend")){
  method<-method[1]
  if(method=="nnls") tree<-nnls.tree(cophenetic(tree),tree,
                                     rooted=TRUE,trace=0)
  else if(method=="extend"){
    h<-diag(vcv(tree))
    d<-max(h)-h
    ii<-sapply(1:Ntip(tree),function(x,y) which(y==x),
               y=tree$edge[,2])
    tree$edge.length[ii]<-tree$edge.length[ii]+d
  } else 
    cat("method not recognized: returning input tree\n\n")
  tree
}

tre.w.names_ultra = force.ultrametric(tre.w.names)
plot(tre.w.names_ultra)



##SSP
SSP_df<- data.frame(names= names(SSP), SSP, groups)
cdat_SSP<- comparative.data(data=SSP_df, phy=tre.w.names_ultra, names.col="names", vcv=TRUE, vcv.dim=3)
#pgls from caper
SSP_mod<- pgls(formula = SSP ~ groups, data = cdat_SSP, lambda='ML')
summary(SSP_mod)
#NS: F-statistic: 2.659 on 1 and 44 DF,  p-value: 0.1101
#AIC(SSP_mod)
resid_ssp<- data.frame(resid(SSP_mod))
resid_ssp_phylo<- phylo4d(tre.w.names, resid_ssp$resid.SSP_mod.)
phyloSignal(resid_ssp_phylo, methods = "K", reps = 10000)



##SSSP
SSSP_df<- data.frame(names= names(SSSP), SSSP, groups)
cdat_SSSP<- comparative.data(data=SSSP_df, phy=tre.w.names_ultra, names.col="names", vcv=TRUE, vcv.dim=3)
#pgls from caper
SSSP_mod<- pgls(formula = SSSP ~ groups, data = cdat_SSSP, lambda='ML')
#SSSP_mod_lambda_kapa_delta  is lowest
summary(SSSP_mod)
#F-statistic: 0.03894 on 1 and 44 DF,  p-value: 0.8445 
resid_sssp<- data.frame(resid(SSSP_mod))
resid_sssp_phylo<- phylo4d(tre.w.names, resid_sssp$resid.SSSP_mod.)
phyloSignal(resid_sssp_phylo, methods = "K", reps = 10000)


##percent SSSP of SSP
SSSP_percent_df<- data.frame(names= names(percentSSSPoutofSSP), percentSSSPoutofSSP, groups)
cdat_SSSP_percent<- comparative.data(data=SSSP_percent_df, phy=tre.w.names_ultra, names.col="names", vcv=TRUE, vcv.dim=3)
#pgls from caper
SSSP_percent_mod<- pgls(formula = percentSSSPoutofSSP ~ groups, data = cdat_SSSP_percent, lambda='ML')
#Choose model
summary(SSSP_percent_mod)
#NS: F-statistic: 0.829 on 1 and 44 DF,  p-value: 0.3675 

#SMC_all
SMC_all_df<- data.frame(names= names(SMC_all), SMC_all, groups)
cdat_SMC_all<- comparative.data(data=SMC_all_df, phy=tre.w.names_ultra, names.col="names", vcv=TRUE, vcv.dim=3)
#pgls from caper
SMC_all_mod<- pgls(formula = SMC_all ~ groups, data = cdat_SMC_all, lambda='ML')
summary(SMC_all_mod)
#NS: F-statistic: 2.467 on 1 and 44 DF,  p-value: 0.1234 

# calculate K of residuals
resid_SMC<- data.frame(resid(SMC_all_mod))
resid_SMC_phylo<- phylo4d(tre.w.names, resid_SMC$resid.SMC_all_mod.)
phyloSignal(resid_SMC_phylo, methods = "K", reps = 10000)


#terp
SMC_terp_df<- data.frame(names= names(SMC_terp), SMC_terp, groups)
cdat_terp<- comparative.data(data=SMC_terp_df, phy=tre.w.names_ultra, names.col="names", vcv=TRUE, vcv.dim=3)
#pgls from caper
terp_mod <- pgls(formula = SMC_terp ~ groups, data = cdat_terp,  lambda='ML')
summary(terp_mod)
#significant: F-statistic: 6.145 on 1 and 44 DF,  p-value: 0.01708 
# calculate K of residuals
resid_terp<- data.frame(resid(terp_mod))
resid_terp_phylo<- phylo4d(tre.w.names, resid_terp$resid.terp_mod.)
phyloSignal(resid_terp_phylo, methods = "K", reps = 10000)



#nrps_like
nrps_df<- data.frame(names= names(SMC_NRPS_like), SMC_NRPS_like, groups)
cdat_nrps<- comparative.data(data=nrps_df, phy=tre.w.names_ultra, names.col="names", vcv=TRUE, vcv.dim=3)
#pgls from caper
nrps_mod <- pgls(formula = SMC_NRPS_like ~ groups, data = cdat_nrps,  lambda='ML')
summary(nrps_mod)
#significant: F-statistic: 5.955 on 1 and 44 DF,  p-value: 0.01878 
resid_nrps<- data.frame(resid(nrps_mod))
resid_nrps_phylo<- phylo4d(tre.w.names, resid_nrps$resid.nrps_mod.)
phyloSignal(resid_nrps_phylo, methods = "K", reps = 10000)


#gpcr
gpcr_df<- data.frame(names= names(GPCRs), GPCRs, groups)
cdat_gpcr<- comparative.data(data=gpcr_df, phy=tre.w.names_ultra, names.col="names", vcv=TRUE, vcv.dim=3)
#pgls from caper
gpcr_mod <- pgls(formula = GPCRs ~ groups, data = cdat_gpcr,  lambda='ML')
summary(gpcr_mod)
#not significant: F-statistic: 0.1227 on 1 and 44 DF,  p-value: 0.7278 
resid_gpcr<- data.frame(resid(gpcr_mod))
resid_gpcr_phylo<- phylo4d(tre.w.names, resid_gpcr$resid.gpcr_mod.)
phyloSignal(resid_gpcr_phylo, methods = "K", reps = 10000)


#get std for table
#all SMCs
std <- function(x) sd(x)/sqrt(length(x))
SMC_traits_S<- SMC_traits[SMC_traits$treatment == "S",]
std_S_SMCs<- std(SMC_traits_S$total_BSG_clusters)
SMC_traits_O<- SMC_traits[SMC_traits$treatment == "O",]
std_S_SMCs<- std(SMC_traits_O$total_BSG_clusters)
  
SMC_traits_red<- SMC_traits_S[SMC_traits_S$host == "R",]
std_SMC_traits_red<- std(SMC_traits_red$total_BSG_clusters)

SMC_traits_white<- SMC_traits_S[SMC_traits_S$host == "W",]
std_SMC_traits_white<- std(SMC_traits_white$total_BSG_clusters)

SMC_traits_larch<- SMC_traits_S[SMC_traits_S$host == "L",]
std_SMC_traits_larch<- std(SMC_traits_larch$total_BSG_clusters)

#terpenes
SMC_traits_S<- SMC_traits[SMC_traits$treatment == "S",]
std_S_SMCs<- std(SMC_traits_S$Terpene)
mean(SMC_traits_S$Terpene)
SMC_traits_O<- SMC_traits[SMC_traits$treatment == "O",]
std_O_SMCs<- std(SMC_traits_O$Terpene)
mean(SMC_traits_O$Terpene)

SMC_traits_red<- SMC_traits_S[SMC_traits_S$host == "R",]
std_SMC_traits_red<- std(SMC_traits_red$Terpene)
mean(SMC_traits_red$Terpene)

SMC_traits_white<- SMC_traits_S[SMC_traits_S$host == "W",]
std_SMC_traits_white<- std(SMC_traits_white$Terpene)
mean(SMC_traits_white$Terpene)

SMC_traits_larch<- SMC_traits_S[SMC_traits_S$host == "L",]
std_SMC_traits_larch<- std(SMC_traits_larch$Terpene)
mean(SMC_traits_larch$Terpene)
std <- function(x) sd(x)/sqrt(length(x))

#NRPS.likes
SMC_traits_S<- SMC_traits[SMC_traits$treatment == "S",]
std_S_SMCs<- std(SMC_traits_S$NRPS.like)
mean(SMC_traits_S$NRPS.like)
SMC_traits_O<- SMC_traits[SMC_traits$treatment == "O",]
std_O_SMCs<- std(SMC_traits_O$NRPS.like)
mean(SMC_traits_O$NRPS.like)

SMC_traits_red<- SMC_traits_S[SMC_traits_S$host == "R",]
std_SMC_traits_red<- std(SMC_traits_red$NRPS.like)
mean(SMC_traits_red$NRPS.like)

SMC_traits_white<- SMC_traits_S[SMC_traits_S$host == "W",]
std_SMC_traits_white<- std(SMC_traits_white$NRPS.like)
mean(SMC_traits_white$NRPS.like)

SMC_traits_larch<- SMC_traits_S[SMC_traits_S$host == "L",]
std_SMC_traits_larch<- std(SMC_traits_larch$NRPS.like)
mean(SMC_traits_larch$NRPS.like)


#get stats without "High" and "Moderate" in Other ECM
#need new tree without the given species 
#high:"Gyrli1" "Rhives1" "Rhivi1"
#moderate: "Rhivul1", Rhisa1
suspects<- c("Gyrli1", "Rhivul1", "Rhivi1", "Rhives1", "Rhisa1")
#shrink tree
tre.w.names_ultra_no_suspects<- drop.tip(phy = tre.w.names_ultra, tip = suspects, trim.internal = TRUE, subtree = FALSE,
         root.edge = 0, rooted = is.rooted(tre.w.names_ultra), collapse.singles = TRUE,
         interactive = FALSE)


##SSP
#remove suspects
SSP_df_wo_suspects<- SSP_df[!rownames(SSP_df) %in% suspects,]
cdat_SSP<- comparative.data(data=SSP_df_wo_suspects, phy=tre.w.names_ultra_no_suspects, names.col="names", vcv=TRUE, vcv.dim=3)
#pgls from caper
SSP_mod <- pgls(formula = SSP ~ groups, data = cdat_SSP, lambda='ML')
summary(SSP_mod)
#not significant

##SSSP
#remove suspects
SSSP_df_wo_suspects<- SSSP_df[!rownames(SSSP_df) %in% suspects,]
cdat_SSSP<- comparative.data(data=SSSP_df_wo_suspects, phy=tre.w.names_ultra_no_suspects, names.col="names", vcv=TRUE, vcv.dim=3)
#pgls from caper
SSSP_mod <- pgls(formula = SSSP ~ groups, data = cdat_SSSP, lambda='ML')
summary(SSSP_mod)
#not significant 

##SMC_all
#remove suspects
SMC_all_df_wo_suspects<- SMC_all_df[!rownames(SMC_all_df) %in% suspects,]
cdat_SMC_all<- comparative.data(data=SMC_all_df_wo_suspects, phy=tre.w.names_ultra_no_suspects, names.col="names", vcv=TRUE, vcv.dim=3)
#pgls from caper
SMC_all_mod <- pgls(formula = SMC_all ~ groups, data = cdat_SMC_all, lambda='ML')
summary(SMC_all_mod)
#not significant overall

##SMC_terp
#remove suspects
SMC_terp_df_wo_suspects<- SMC_terp_df[!rownames(SMC_terp_df) %in% suspects,]
cdat_SMC_terp<- comparative.data(data=SMC_terp_df_wo_suspects, phy=tre.w.names_ultra_no_suspects, names.col="names", vcv=TRUE, vcv.dim=3)
#pgls from caper
SMC_terp_mod <- pgls(formula = SMC_terp ~ groups, data = cdat_SMC_terp, lambda='ML')
summary(SMC_terp_mod)
#Residual standard error: 9.692 on 39 degrees of freedom
#Multiple R-squared: 0.4275,	Adjusted R-squared: 0.4128 
#F-statistic: 29.12 on 1 and 39 DF,  p-value: 3.567e-06 

##SMC_nrps
#remove suspects
SMC_nrps_df_wo_suspects<- nrps_df[!rownames(nrps_df) %in% suspects,]
cdat_SMC_nrps<- comparative.data(data=SMC_nrps_df_wo_suspects, phy=tre.w.names_ultra_no_suspects, names.col="names", vcv=TRUE, vcv.dim=3)
#pgls from caper
SMC_nrps_mod <- pgls(formula = SMC_NRPS_like ~ groups, data = cdat_SMC_nrps, lambda='ML')
summary(SMC_nrps_mod)
#Residual standard error: 5.149 on 39 degrees of freedom
#Multiple R-squared: 0.2241,	Adjusted R-squared: 0.2042 
#F-statistic: 11.27 on 1 and 39 DF,  p-value: 0.001771 


##gpcr
#remove suspects
gpcr_df_wo_suspects<- gpcr_df[!rownames(gpcr_df) %in% suspects,]
cdat_gpcr<- comparative.data(data=gpcr_df_wo_suspects, phy=tre.w.names_ultra_no_suspects, names.col="names", vcv=TRUE, vcv.dim=3)
#pgls from caper
gpcr_mod <- pgls(formula = GPCRs ~ groups, data = cdat_gpcr, lambda='ML')
summary(gpcr_mod)
#not significant

