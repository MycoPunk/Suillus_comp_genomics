#Running R version 3.5.3

#set directory
setwd("~/Desktop")

#load libs
library("phytools")
library("ape")
library("quadprog")
library("devtools")
library("phangorn")

state_df_stag<- read.csv("Suillus_trait_states_STAG.csv")
state_df_iq<- read.csv("Suillus_trait_states_IQ.csv")

#read in tree
phy_stag<- as.character(state_df_stag[1,5])
phy_stag<- read.tree(text = phy_stag)

phy_iq<- as.character(state_df_iq[1,5])
phy_iq<- read.tree(text = phy_iq)

#read in the state data
host_state_df <- data.frame(taxa=state_df_stag$X, host_state=(state_df_stag[,4]))

x2<- host_state_df[,2] 
x3<-setNames(host_state_df[,2],host_state_df[,1])

#flip all
phy_stag<-rotateNodes(phy_stag,"all")
phy_iq<-rotateNodes(phy_iq,"all")

#need to make a simmap to connect chr states to tips
phy_stag<-make.simmap(phy_stag, x3)
phy_iq<-make.simmap(phy_iq, x3)

#link the host states
states_stag<-getStates(phy_stag,"tips")
states_iq<-getStates(phy_iq,"tips")

#set color pallet: 
#(#3E4B60 = dark blue = generalist)
#(#FFFFFF) = white = unkown 
#(#B36757 = pink = Larch)
#(#405952 = dark green = Pseudotsuga)
#(#5676A1 = light blue = Red Pine)
#(#B09136 = gold = White Pine)
#palette = c("#3E4B60", "#B36757", "#405952", "#5676A1", "#B09136")

#plot IQ tree
cols<-setNames(palette(c("#3E4B60", "#B36757", "#405952", "#5676A1", "#B09136", "#FFFFFF"))[1:length(unique(states_iq))],sort(unique(states_iq)))
fitER<-ace(states_iq,phy_iq,model="ER",type="discrete")
fitER
plotTree(phy_iq,type="phylogram",fsize=0.8,ftype="i")
nodelabels(node=1:phy_iq$Nnode+Ntip(phy_iq),
           pie=fitER$lik.anc,piecol=cols,cex=0.5)
tiplabels(pie=to.matrix(states_iq,sort(unique(states_iq))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(phy_iq)),fsize=0.8)


#plot STAG tree
cols2<-setNames(palette(c("#3E4B60", "#B36757", "#405952", "#5676A1", "#FFFFFF", "#B09136"))[1:length(unique(states_stag))],sort(unique(states_stag)))
fitER<-ace(states_stag,phy_stag,model="ER",type="discrete")
fitER
plotTree(phy_stag,type="phylogram",fsize=0.8,ftype="i")
nodelabels(node=1:phy_stag$Nnode+Ntip(phy_stag),
           pie=fitER$lik.anc,piecol=cols2,cex=0.5)
tiplabels(pie=to.matrix(states_stag,sort(unique(states_stag))),piecol=cols2,cex=0.3)
add.simmap.legend(colors=cols2,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(phy_stag)),fsize=0.8)


#create object cophylo and rotate nodes to optimize matching 
matched_trees<-cophylo(phy_iq,phy_stag, rotate=TRUE)
plot(matched_trees,fsize=0.5)
#gives error  "Error in cat(x, file = file, sep = c(rep.int(sep, ncolumns - 1), "\n"),  : 
#  invalid connection"


#try resetting all the inputs and creating cophylo object wihout linked trait states
#re-read this all overwrite vars:
#read in tree
phy_stag<- as.character(state_df_stag[1,5])
phy_stag<- read.tree(text = phy_stag)
phy_iq<- as.character(state_df_iq[1,5])
phy_iq<- read.tree(text = phy_iq)


#flip all
phy_stag<-rotateNodes(phy_stag,"all")
phy_iq<-rotateNodes(phy_iq,"all")

matched_trees<-cophylo(phy_iq,phy_stag, rotate=TRUE)
plot(matched_trees,fsize=0.5)

#make a simmap to connect chr states to tips
phy_iq<-make.simmap(phy_iq, states_iq)
phy_stag<-make.simmap(phy_stag, states_stag)

#link the host states
states_iq<-getStates(phy_iq,"tips")
states_stag<-getStates(phy_stag,"tips")



##choose model IQ
# equal rates model 
fitER<-fitMk(phy_stag,states_stag,model="ER")
fitER
plot(fitER)

# SYM
fitSYM<-fitMk(phy_stag,states_stag,model="SYM")
fitSYM
plot(fitSYM, show.zeros=FALSE)

# ARD 
fitARD<-fitMk(phy_stag,states_stag,model="ARD")
fitARD
plot(fitARD,show.zeros=FALSE)

#get AIC
AIC_stag<-setNames(sapply(list(fitER,fitSYM,fitARD),AIC),c("ER","SYM","ARD"))
AIC_stag
aic.w(AIC_stag)
#ER wins STAG
#ER        SYM        ARD 
#0.85389636 0.14605968 0.00004397 

##choose model IQ
# equal rates model 
fitER<-fitMk(phy_iq,states_iq,model="ER")
fitER
plot(fitER)

# SYM
fitSYM<-fitMk(phy_iq,states_iq,model="SYM")
fitSYM
plot(fitSYM, show.zeros=FALSE)

# ARD 
fitARD<-fitMk(phy_iq,states_iq,model="ARD")
fitARD
plot(fitARD,show.zeros=FALSE)

#get AIC
AIC_iq<-setNames(sapply(list(fitER,fitSYM,fitARD),AIC),c("ER","SYM","ARD"))
AIC_iq
aic.w(AIC_iq)
#ER wins Siq
#ER        SYM        ARD 
#0.84556852 0.15433521 0.00009626


#set node colors
cols<-setNames(palette(c("#3E4B60", "#B36757", "#405952", "#5676A1", "#B09136"))[1:length(unique(states_iq))],sort(unique(states_iq)))

#make models
model_iq<-ace(states_iq,matched_trees$trees[[1]],model="ER",type="discrete")
model_stag<-ace(states_stag,matched_trees$trees[[2]],model="ER",type="discrete")
#plot trees and ancesteral states
plot(matched_trees, link.type="curved", fsize=c(0.6,0.6))
nodelabels.cophylo(pie=model_iq$lik.anc,piecol=cols, which="left", cex = 0.6)
nodelabels.cophylo(pie=model_stag$lik.anc,piecol=cols, which="right", cex = 0.6)
tiplabels.cophylo(pie=to.matrix(states_iq[matched_trees$trees[[1]]$tip.label],sort(unique(states_iq))), piecol=cols,cex=0.4, which="left")
tiplabels.cophylo(pie=to.matrix(states_stag[matched_trees$trees[[1]]$tip.label],sort(unique(states_stag))), piecol=cols,cex=0.4, which="right")


