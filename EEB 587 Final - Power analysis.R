# EEB 587 Final Project - due the last day of finals (18 May 2022)

# 1. Create project question 
#Power analysis to decide the appropriate number/range of taxa needed for dissertation work 

# Question 1 - What is an appropriate number of taxa to use when estimating the co-variation of two traits?

# Question 2 - Does this minimum number of taxa differ when using two different primate phylogenies? 

#Methods/steps:
# Test the Boettiget, Coop, and Ralph, 2012 (Data from: Is your phylogeny informative? Measuring the power of comparative methods) package 

#2. Input and read primate tree 
library(ape)
library(geiger)
library(phytools)
library(Runuran)

primate.tree.test <- "TreeBlock_10kTrees_Primates_Version3.nex"
ape::read.nexus(primate.tree.test)

primate.MF.tree <- "primate_trees.nex copy"
ape::read.nexus(primate.MF.tree)
lapply(ape::read.nexus(primate.MF.tree, force.multi = TRUE), ape::consensus)

#Use the following tree:
#Downloaded from the 10k Trees Primates program 
#Branch lengths are proportional to absolute time. Ultrametric tree as Chronogram 
primate.tree2 <- ape::read.nexus("consensusTree_10kTrees_Primates_Version3.nex")
plot(primate.tree2)
primate.tree2$edge
str(primate.tree2) #structure of the tree
primate.tree2$tip.label #vector of names 



# 3. Simulate data
# Want to estimate co-variation between two traits 
# Possible co-variation values:
possible_covar <- c(0, 0.5, 1)

#Trait 1 (Long-axis length of thyrohyal [greater horn] of hyoid) 10 mm
Trait_1 <- 10

#Trait 2 (Length of clavicle) 145 mm 
Trait_2 <- 145

# Generate random data:
set.seed(1234)
x <- rnorm(n=301, mean=10, sd=5)
require("Runuran")
data.T1 <- urnorm(n=301, mean=10, sd=2, ub=100)
data.T2 <- urnorm(n=301, mean= 145, sd=25, ub=100)
dt <- cbind(data.T1, data.T2)
data.primate <- cbind(primate.tree2$tip.label, data.T1, data.T2)
print(as.data.frame(sapply(data.primate, function(x) gsub("\"","",x))))
primate.data.edited <- data.frame(primate.tree2$tip.label, dt, stringsAsFactors = FALSE); print(primate.data.edited, quote=FALSE)


# 4. Figure out a good sigma^2 
test1_sigma <- c(0, 0.1, 0.2)

#Number of taxa to keep 
possible_taxa <- round(Ntip(primate.tree2)*c(.1, .5, 1))
str(possible_taxa)

covar_index <- seq(length(possible_covar))
sigma_index <- seq(length(test1_sigma))
replicate_index <- sequence(10)
covar_matrix <- matrix(test1_sigma[sigma_index], nrow=3, ncol=3)
covar_matrix[1,2] <- possible_covar[covar_index]
covar_matrix[2,1] <- covar_matrix[1,2]
traits <- sim.char(primate.tree2, par = covar_matrix, nsim = 1, model = "BM", root=c(Trait_1, Trait_2))
sim.traits <- sim.char(primate.tree2, 0.02, 100)
sim.trait.test <- sim.char(primate.tree2, 0.2, 100)

taxa_index <- sequence(possible_taxa)
phy_pruned <- geiger::drop.random(primate.tree2, n=ape::Ntip(primate.tree2))
str(phy_pruned)
pruned_all <- geiger::treedata(phy_pruned, sim.trait.test)
geiger_results <- geiger::fitContinuous(phy_pruned$primate.tree2, dat=primate.data.edited)



# 5. Figure out the likelihood 
BM.fit.lik <- fitContinuous(phy=primate.tree2, dat=primate.data.edited, SE=NA, control = list(niter=50), ncores = 2)
bm2 <- fitContinuous(phy=primate.tree2, dat=data.primate, SE=NA, control = list(niter=50), ncores = 2)
# 6. Analyze the simulations 


# 7. Visualize 
library(tidyverse)
library(ggplot2)
library(ggtree)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ggtree")


testprime <- primate.tree2
d1 <- data.frame(id=testprime$tip.label, val=primate.data.edited[,2])
p2 <- facet_plot(testprime, panel= "dot", data=d1, geom=geom_point, aes(x=val), color='red3')
d2 <- data.frame(id=primate.tree2$tip.label, value=primate.data.edited[,3])
p3 <- facet_plot(p2, panel = 'bar', data = d2, geom = geom_segment, aes(x=0, xend=value, y=y, yend=y), size =3, color='blue4')
p3 + theme_tree2()

# 8. Decide the minimum number of taxa needed/appropriate 

# 9. Compare results from two different primate tree 
# Comparision tree from 10k Trees Primate Project 
primate.tree3 <- ape::read.nexus("consensusTree_10kTrees_Primates_Version3_withCladeCredibilityValues.nex")
plot(primate.tree2)












































































