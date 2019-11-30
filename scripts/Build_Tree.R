#set directory
setwd("~/Desktop/Project_Suillus_comp_genomics/R")
library("phytools")
library(ape)
library("quadprog")
install.packages("quadprog")
options(stringsAsFactors = TRUE)
library(devtools)
install_github("liamrevell/phytools")
biocLite("Biostrings")
R.version
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

source("https://bioconductor.org/biocLite.R")

#load tree wih branch lengths generate with STAG
#note- there is a bug in this protram, where pies representing ancesteral state give an error for $ based identification.
#if it does this, run the example annole script at the bottom of the page, and the pies start working again (I dono). 

state_df<- read.csv("Suillus_trait_states_pub.csv")

#read in tree
phy<- as.character(state_df[1,6])
phy<- read.tree(text = phy)

plot.phylo(phy)


#read in the data
host_state_df <- data.frame(taxa=state_df$X, host_state=(state_df[,5]))


x2<- host_state_df[,2] 
x3<-setNames(host_state_df[,2],host_state_df[,1])


#need to make a simmap to connect chr states to tips
phy2<-make.simmap(phy, x3)

#link the host states - this can take a bit of time
states<-getStates(phy2,"tips")
states

#this is your tree
plotTree(phy2,type="phylogram",fsize=0.8,ftype="i")

#track states no coloring the nodes
cols<-setNames(palette()[1:length(unique(states))],sort(unique(states)))
tiplabels(pie=to.matrix(states,sort(unique(states))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(phy2)),fsize=0.8)


#now calculate probability of ancesteral state
fitER<-ace(states,phy2,model="ER",type="discrete")
fitER
plotTree(phy2,type="phylogram",fsize=0.8,ftype="i")
nodelabels(node=1:tree$Nnode+Ntip(phy2),
           pie=fitER$lik.anc,piecol=cols,cex=0.5)
tiplabels(pie=to.matrix(states,sort(unique(states))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(phy2)),fsize=0.8)



#the above but using an MCMC approach to sample character histories from their posterior probability distribution

# simulate single stochastic character map using empirical Bayes method
mtree<-make.simmap(phy2,states,model="ER")
mtree
plot(mtree,type="phylogram",fsize=0.8,ftype="i")
add.simmap.legend(colors=cols,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(phy2)),fsize=0.8)


#Note- this dosent really mean anything, we need to look at a whole population of these mappings
mtrees<-make.simmap(phy2,states,model="ER",nsim=100)
mtrees

par(mfrow=c(10,10))
null<-sapply(mtrees,plot,colors=cols,lwd=1,ftype="off")

#yuck. lets turn that into stats
#this is the proportion of time spent in each state, and the posterior probabilities that each internal node is in each state, under the model.
pd<-summary(mtrees,plot=FALSE)
pd


par(mfrow = c(1,1)) #back to one image per page


# #example data for other trees you can build in this package

data(anoletree)
## this is just to pull out the tip states from the
## data object - normally we would read this from file
x<-getStates(anoletree,"tips")
tree<-anoletree
rm(anoletree)
tree
x
plotTree(tree,type="fan",fsize=0.8,ftype="i")
cols<-setNames(palette()[1:length(unique(x))],sort(unique(x)))
tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(tree)),fsize=0.8)

fitER<-ace(x,tree,model="ER",type="discrete")
fitER
plotTree(tree,type="fan",fsize=0.8,ftype="i")
nodelabels(node=1:tree$Nnode+Ntip(tree),
           pie=fitER$lik.anc,piecol=cols,cex=0.5)
tiplabels(pie=to.matrix(x,sort(unique(x))),piecol=cols,cex=0.3)
add.simmap.legend(colors=cols,prompt=FALSE,x=0.9*par()$usr[1],
                  y=-max(nodeHeights(tree)),fsize=0.8)


mtree<-make.simmap(tree,x,model="ER")
mtree
plot(mtree,cols,type="fan",fsize=0.8,ftype="i")