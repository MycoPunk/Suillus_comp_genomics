#this script investigates CAZy enzymes, to compare degredation ability in Suillus vs. other ECM fungi. 
#with special interest in the role of Auxiliary Activitie groups, putative involved in SOM oxidation.

#set wd
setwd("~/Desktop/")

#set libraries
library(phylosignal)
library(adephylo)
library(ape)
library(phylobase)
library(phytools)
library(caper)
library(geiger)
library(pheatmap)
library(tidyr)
library(dplyr)
library(phangorn) 
library(tidyverse)

#set seed for reproducibility
set.seed(666)

##read in trait data
CAZy_traits<- read.table(file ="CAZy_data.csv", header = TRUE, sep = ",")

#make the .phylo file
#read in the tree
tree<- "((((((((((((((((((Suivar1:0.010372,Suitom1:0.006211):0.001161,(Suihi1:0.009804,Suifus1:0.009871):0.001006):0.003248,Suiplo1:0.010465):0.00292,Suidis1:0.010283):0.018954,Suibov1:0.042762):0.011931,(((((Suisubl1:0.010257,Suisu1:0.015162):0.013213,Suiame1:0.020213):0.003498,Suisub1:0.0274):0.003311,Suicot1:0.032812):0.00549,(Suipic1:0.020678,Suidec1:0.01842):0.024794):0.003182):0.0023,(Suigr1:0.01368,Suipla1:0.010526):0.016982):0.002008,((Suilu4:0.010432,Suiocc1:0.016253):0.004298,Suibr2:0.011667):0.022327):0.010157,Suilak1:0.033981):0.008243,Suicli1:0.034744):0.017466,(Suipal1:0.037565,Suiamp1:0.034194):0.010046):0.025448,((((Rhivi1:0.013169,Rhives1:0.018102):0.028885,Rhivul1:0.057818):0.006656,Rhisa1:0.051235):0.007496,Rhitru1:0.079972):0.022267):0.113665,(((Pismi1:0.067099,Pisti1:0.055341):0.081864,Sclci1:0.12731):0.085541,(Paxin1:0.06372,Gyrli1:0.062439):0.081453):0.055982):0.104422,Pilcr1:0.205084):0.029249,(((Lacam2:0.030668,Lacbi2:0.031522):0.145046,Hebcy2:0.199426):0.059435,Amamu1:0.281571):0.070896):0.042081,((Theter1:0.080412,Thega1:0.073159):0.28376,(Ruscom1:0.087985,Rusbre1:0.092679):0.251688):0.035305):0.111335,(Hydru2:0.235556,Cananz1:0.352162):0.263901):0.092299,(Gaumor1:0.166361,Hyssto1:0.187679):0.092299);"
tre.w.names<- read.tree(text = tree)
#make phylogenetic distance matrix
tre.w.names_coph<- cophenetic(tre.w.names)

#replace NAs with zeros
CAZy_traits[is.na(CAZy_traits)] <- 0


#are any of the categores significantly different between the groups? 
#isolate trait data primary categories (AA, CBM, CE, EXPN, GH162, GH, GT, Myosin motor,PL)
AA<- CAZy_traits$AA
names(AA) = CAZy_traits$JGI_project_code

CBM<- CAZy_traits$CBM
names(CBM) = CAZy_traits$JGI_project_code

CE<- CAZy_traits$CE
names(CE) = CAZy_traits$JGI_project_code

EXPN<- CAZy_traits$EXPN
names(EXPN) = CAZy_traits$JGI_project_code

GH162<- CAZy_traits$GH162
names(GH162) = CAZy_traits$JGI_project_code

GH<- CAZy_traits$GH
names(GH) = CAZy_traits$JGI_project_code

GT<- CAZy_traits$GT
names(GT) = CAZy_traits$JGI_project_code

Myosin_motor<- CAZy_traits$Myosin_motor
names(Myosin_motor) = CAZy_traits$JGI_project_code

PL<- CAZy_traits$PL
names(PL) = CAZy_traits$JGI_project_code


#order the trait values to match the tree order AA, CBM, CE, EXPN, GH162, GH, GT, Myosin motor,PL
AA<- AA[row.names(tre.w.names_coph)]
CBM<- CBM[row.names(tre.w.names_coph)]
CE<- CE[row.names(tre.w.names_coph)]
EXPN<- EXPN[row.names(tre.w.names_coph)]
GH162<- GH162[row.names(tre.w.names_coph)]
GH<- GH[row.names(tre.w.names_coph)]
Myosin_motor<- Myosin_motor[row.names(tre.w.names_coph)]
PL<- PL[row.names(tre.w.names_coph)]

#set gorups
groups<- CAZy_traits$treatment
names(groups) = CAZy_traits$JGI_project_code
#order the trait values to match the tree order
groups<- groups[row.names(tre.w.names_coph)]



###phylosgignal graphs for local autocorrelation
#AA, CBM, CE, EXPN, GH162, GH, GT, Myosin motor,PL
CAZy_input<- data.frame(AA = AA, CBM = CBM, CE = CE, EXPN = EXPN, GH162 = GH162, GH = GH, GT = GT, Myosin_motor = Myosin_motor, PL = PL)
#create .p4d object
CAZy_input.p4d_named <- phylo4d(tre.w.names, CAZy_input)
#compute lipaMoran I 
CAZy.lipa <- lipaMoran(CAZy_input.p4d_named, reps = 10000)
CAZy.lipa.p4d <- lipaMoran(CAZy_input.p4d_named, as.p4d = TRUE, reps = 10000)
#set colors and cut offs
points.col<-ifelse(CAZy.lipa$p.value < 0.01, "#FCD090", "#93959799")
#plot -AA
barplot.phylo4d(CAZy_input.p4d_named, bar.col=points.col, 
                center = TRUE, 
                scale = FALSE, 
                tree.ladderize = TRUE,
                tip.cex = .5,
                trait.cex = .5,
                trait.bg.col = "white",
                grid.vertical = FALSE,
                grid.horizontal = FALSE,
                bar.lwd = 5,
                rotate.tree = 90)

#make tree ultimetric for pgls analysis 
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


##PGLS on all 
#step 1 sort the main traits df
CAZy_traits_ordered<-CAZy_traits[match(row.names(tre.w.names_coph), CAZy_traits$JGI_project_code),]

#sort groups too 
groups_ordered<- groups[row.names(tre.w.names_coph)]

#format the dataframe
CAZy_traits2<-as.matrix(CAZy_traits_ordered[,4:ncol(CAZy_traits_ordered)])
row.names(CAZy_traits2)<- CAZy_traits_ordered$JGI_project_code
CAZy_traits3<- data.frame(CAZy_traits2)
#transform
CAZy_traits3_t<- t(CAZy_traits3)
#split the df by row into a list of lists
data_list <- setNames(split(CAZy_traits3_t,
                            seq(nrow(CAZy_traits3_t))),
                      rownames(CAZy_traits3_t))
#name each element in the list
data_list_named<- lapply(data_list, function(x) setNames(x, CAZy_traits_ordered$JGI_project_code))


#step 2: make each column into a dataframe
#form this: AA_df<- data.frame(names= names(AA), AA, groups)
CAZY_df<- data.frame()
CAZY_df_list<- vector(mode = "list", length = length(data_list_named))
for (j in 1:length(data_list_named)){
  for (i in data_list_named[j]){
    CAZY_df <-  data.frame(names= names(i), i, groups)
    CAZY_df_list[[j]]<- CAZY_df
  }
}

#name ea. df in the list of dfs
names(CAZY_df_list) = names(data_list_named)


##step 2: link each dataframe to the phylogeny creating a "comparative dataset"
# from this: cdat_AA<- comparative.data(data=AA_df, phy=tre.w.names_ultra, names.col="names", vcv=TRUE, vcv.dim=3)

#note- this works but takes a hot minute 
list_of_cdat<- vector(mode = "list", length = length(data_list_named))
for (j in 1:length(list_of_cdat)){
  for (i in CAZY_df_list[j]){
    ea_cdat<- comparative.data(data=i, phy=tre.w.names_ultra, names.col="names", vcv=TRUE, vcv.dim=3)
    list_of_cdat[[j]]<- ea_cdat
  }}

#name ea. df in the list of dfs
names(list_of_cdat) = names(data_list_named)


#step 3: run pgls
#pgls from caper
#from this: AA_mod<- pgls(formula = AA ~ groups, data = cdat_AA, lambda='ML')
#summary(AA_mod)
list_of_mod<- vector(mode = "list", length = length(data_list_named))
for (j in 1:length(list_of_mod)){
  for (i in list_of_cdat[j]){
    mod_ea<- summary(pgls(formula = i ~ groups, data = i, lambda='ML'))
    list_of_mod[[j]]<- mod_ea
  }}

#name ea. df in the list of dfs
names(list_of_mod) = names(data_list_named)


#extract p-values 
list_of_p_vals<- vector(mode = "list", length = length(data_list_named))
for (j in 1:length(list_of_p_vals)){
  for (i in list_of_mod[j]){
    pval_ea<- i$coefficients[2,4]
    list_of_p_vals[[j]]<- pval_ea
  }}
#name ea. df in the list of dfs
names(list_of_p_vals) = names(data_list_named)

#how many are significant?
length(list_of_p_vals[list_of_p_vals <.01])
#there are 28 that are significantly different
#there are 12 that are significantly different at p <0.01

#the different ones:
different_between_SandO<- list_of_p_vals[list_of_p_vals <.01]



##remove suspects and re-run to look for changes in significant groups
suspects<- c("Gyrli1", "Rhivul1", "Rhivi1", "Rhives1", "Rhisa1")
CAZy_traits_ordered_no_suspects<- CAZy_traits_ordered[!CAZy_traits_ordered$JGI_project_code %in% suspects, ]
#remove suspects from groups
groups_no_suspects<- groups[!names(groups) %in% suspects]
#remove suspects from tree
tre.w.names_no_suspects<- drop.tip(tre.w.names, suspects, trim.internal = TRUE, subtree = FALSE,
                              root.edge = 0, collapse.singles = TRUE,
                              interactive = FALSE)
#make cophenetic matrix
tre.w.names_no_suspects_coph<- cophenetic(tre.w.names_no_suspects)
#make tree ultametric
tre.w.names_ultra_ns = force.ultrametric(tre.w.names_no_suspects)

#check that the order of the dataframe matches the tree now that you've removed some tips
CAZy_traits_ordered_no_suspects_again<-CAZy_traits_ordered_no_suspects[match(row.names(tre.w.names_no_suspects_coph), CAZy_traits_ordered_no_suspects$JGI_project_code),]

#format the dataframe
CAZy_traits_ordered_no_suspects2<-as.matrix(CAZy_traits_ordered_no_suspects_again[,4:ncol(CAZy_traits_ordered_no_suspects_again)])
row.names(CAZy_traits_ordered_no_suspects2)<- CAZy_traits_ordered_no_suspects_again$JGI_project_code
CAZy_traits_ordered_no_suspects3<- data.frame(CAZy_traits_ordered_no_suspects2)

#transform
CAZy_traits_ordered_no_suspects3_t<- t(CAZy_traits_ordered_no_suspects3)
#split the df by row into a list of lists
data_list_ns <- setNames(split(CAZy_traits_ordered_no_suspects3_t,
                            seq(nrow(CAZy_traits_ordered_no_suspects3_t))),
                      rownames(CAZy_traits_ordered_no_suspects3_t))
#name each element in the list
data_list_named_ns<- lapply(data_list_ns, function(x) setNames(x, CAZy_traits_ordered_no_suspects_again$JGI_project_code))


#step 2: make each column into a dataframe
#form this: AA_df<- data.frame(names= names(AA), AA, groups)
CAZY_df_ns<- data.frame()
CAZY_df_list_ns<- vector(mode = "list", length = length(data_list_named_ns))
for (j in 1:length(data_list_named_ns)){
  for (i in data_list_named_ns[j]){
    CAZY_df_ns <-  data.frame(names= names(i), i, groups_no_suspects)
    CAZY_df_list_ns[[j]]<- CAZY_df_ns
  }
}

#name ea. df in the list of dfs
names(CAZY_df_list_ns) = names(data_list_named_ns)


##step 2: link each dataframe to the phylogeny creating a "comparative dataset"
# from this: cdat_AA<- comparative.data(data=AA_df, phy=tre.w.names_ultra, names.col="names", vcv=TRUE, vcv.dim=3)

#note- this works but takes a hot minute 
list_of_cdat_ns<- vector(mode = "list", length = length(data_list_named_ns))
for (j in 1:length(list_of_cdat_ns)){
  for (i in CAZY_df_list_ns[j]){
    ea_cdat_ns<- comparative.data(data=i, phy=tre.w.names_ultra_ns, names.col="names", vcv=TRUE, vcv.dim=3)
    list_of_cdat_ns[[j]]<- ea_cdat_ns
  }}

#name ea. df in the list of dfs
names(list_of_cdat_ns) = names(data_list_named_ns)



#step 3: run pgls
#pgls from caper
#from this: AA_mod<- pgls(formula = AA ~ groups, data = cdat_AA, lambda='ML')
#summary(AA_mod)
list_of_mod_ns<- vector(mode = "list", length = length(data_list_named_ns))
for (j in 1:length(list_of_mod_ns)){
  for (i in list_of_cdat_ns[j]){
    mod_ea_ns<- pgls(formula = i ~ groups_no_suspects, data = i, lambda='ML')
    list_of_mod_ns[[j]]<- mod_ea_ns
  }}
#note - this throws a weird error: Error in pgls(formula = i ~ groups_no_suspects, data = i, lambda = "ML") : 
#Problem with optim:52ERROR: ABNORMAL_TERMINATION_IN_LNSRCH
#go back and do this individually if needed - it's probably due to 0 variance in some comparisons 
#if you're only looking at AA's this isn't hard to do and will reduce the chance to running into this model error

#name ea. df in the list of dfs
names(list_of_mod_ns) = names(data_list_named_ns)


#extract p-values 
list_of_p_vals<- vector(mode = "list", length = length(data_list_named))
for (j in 1:length(list_of_p_vals)){
  for (i in list_of_mod[j]){
    pval_ea<- i$coefficients[2,4]
    list_of_p_vals[[j]]<- pval_ea
  }}
#name ea. df in the list of dfs
names(list_of_p_vals) = names(data_list_named)

#how many are significant?
length(list_of_p_vals[list_of_p_vals <.05])
#there are 28 that are significantly different

#the different ones:
different_between_SandO<- list_of_p_vals[list_of_p_vals <.05]

##make a heatmap for the CAZy's that are  significantly different 
sig_diff_names<- names(different_between_SandO) 
#+ remove the overall GH category 
sig_diff_names2<-sig_diff_names[unlist(sig_diff_names) != "GH"]
#subset the input to just these results:
CAZy_traits3_significant<- CAZy_traits3[,colnames(CAZy_traits3) %in% sig_diff_names2]

#slide that over to GHs
p_sig <-  gheatmap(p8, CAZy_traits3_significant, 
                 offset=0.3, 
                 width=2, 
                 low="white", 
                 high="black", 
                 colnames_position = "top", 
                 font.size=1.5,
                 colnames_angle=45,
                 colnames_offset_y = 1)

# plot
plot(p_sig)

#in color useing pheatmap
#order the dataframe 
names_in_order<- c("Suibov1",
                   "Suibr2",
                   "Suicot1",
                   "Suidec1",
                   "Suifus1",
                   "Suihi1",
                   "Suilu4",
                   "Suiocc1",
                   "Suisu1",
                   "Suitom1",
                   "Suivar1",
                   "Suiame1",
                   "Suidis1",
                   "Suipic1",
                   "Suipla1",
                   "Suiplo1",
                   "Suisubl1",
                   "Suiamp1",
                   "Suicli1",
                   "Suipal1",
                   "Suisub1",
                   "Suigr1",
                   "Suilak1",
                   "Rhisa1",
                   "Rhitru1",
                   "Rhives1",
                   "Rhivi1",
                   "Rhivul1",
                   "Amamu1",
                   "Cananz1",
                   "Gaumor1",
                   "Gyrli1",
                   "Hebcy2",
                   "Hydru2",
                   "Hyssto1",
                   "Lacam2",
                   "Lacbi2",
                   "Paxin1",
                   "Pilcr1",
                   "Pismi1",
                   "Pisti1",
                   "Rusbre1",
                   "Ruscom1",
                   "Sclci1",
                   "Theter1",
                   "Thega1")


CAZy_traits3_significant_ordered<- CAZy_traits3_significant[match(names_in_order, rownames(CAZy_traits3_significant)),]   

pheatmap(CAZy_traits3_significant_ordered, 
         angle_col = "315", 
         fontsize = 3,
         cluster_rows= FALSE,
         cluster_cols = TRUE,
         scale = "column",
         color = inferno(100),
         #color = colfunc(100),
         #breaks = 1)
         border_color = NA,
         cellheight=5,cellwidth=5,
         display_numbers = CAZy_traits3_significant_ordered)

CAZy_traits3_AA_ordered<- CAZy_traits3[match(names_in_order, rownames(CAZy_traits3)),]   

#do the same but with all AA 
pheatmap(CAZy_traits3_AA_ordered[2:21], 
         angle_col = "315", 
         fontsize = 3,
         cluster_rows= FALSE,
         cluster_cols = TRUE,
         scale = "column",
         color = inferno(100),
         #color = colfunc(100),
         #breaks = 1)
         border_color = NA,
         cellheight=5,cellwidth=5,
         display_numbers = CAZy_traits3_AA_ordered[2:21])

#get average 
AA3_dist_df<- data.frame(cbind(CAZy_traits$JGI_project_code, CAZy_traits$treatment, as.numeric(CAZy_traits$AA3_dist)))
sul<- AA3_dist_df[AA3_dist_df$X2 == "S",]
oth<- AA3_dist_df[AA3_dist_df$X2 == "O",]
mean(as.numeric(sul$X3))
mean(as.numeric(oth$X3))

#get pval for AA3_dist
list_of_mod$AA3_dist

