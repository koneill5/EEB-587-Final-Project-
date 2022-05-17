# EEB 587 Final Project - due the last day of finals (18 May 2022)

# 1. Create project question 
#Power analysis to decide the appropriate number/range of taxa needed for dissertation work 

# Question 1 - What is an appropriate number of taxa to use when estimating the co-variation of two traits?

# Question 2 - Does this minimum number of taxa differ when using two different primate phylogenies? 

#Methods/steps:

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

#### Next Section ####

# 3. Simulate data
# Want to estimate co-variation between two traits 
# Possible co-variation values:
possible_covar <- c(0, 0.5, 1)

#Trait 1 (Long-axis length of thyrohyal [greater horn] of hyoid) average: 10 mm
Trait_1 <- 10

#Trait 2 (Length of clavicle) average: 145 mm 
Trait_2 <- 145

# Generate random data:
set.seed(1234)
x <- rnorm(n=301, mean=10, sd=5)
require("Runuran")
data.T1 <- urnorm(n=301, mean=10, sd=2,lb=1, ub=100)
data.T2 <- urnorm(n=301, mean= 145, sd=25,lb=1, ub=200)
dt <- cbind(data.T1, data.T2)
data.primate <- cbind(primate.tree2$tip.label, data.T1, data.T2)
data.primate1 <- as.data.frame(data.primate)
name.check(primate.tree2,data.primate1)

### Trying to fix the quotation issue: 
print(as.data.frame(sapply(data.primate1, function(x) gsub("\"","",x))))
primate.data.edited <- data.frame(primate.tree2$tip.label, data.primate1, stringsAsFactors = TRUE); print(primate.data.edited, quote=FALSE)
primate.data.edited$primate.tree2.tip.label<-paste0("",primate.data.edited$primate.tree2.tip.label,"")

name.check(primate.tree2,primate.data.edited$primate.tree2.tip.label)

##### Next section ####

# Compare taxa in dat and tree:
treedata(primate.tree2, primate.data.edited, sort=FALSE, warnings=TRUE)
tmp <- name.check(primate.tree2, primate.data.edited, data.names = NULL)
newphy <- drop.tip(primate.tree2, tip=tmp$tree_not_data)
name.check(newphy,primate.data.edited)

# Match names:
new.phy <- primate.data.edited[,1]
new.phy2 <- sub("", "", new.phy)
names(primate.tree2.tip.label) <- new.phy2


#### Next section ####

#Map continuous traits on the tree 
#Extract character of interest: 
horn.length <- setNames(primate.data.edited$data.T1, rownames(primate.data.edited[,1]))

# brownie lite
multiBM.fit.test <- brownie.lite(primate.tree2,horn.length)


#Contmap 
primate.contmap <- contMap(primate.tree2, horn.length, plot = FALSE, res = 200)

#### Next section ####

# 4. Figure out a good sigma^2 
test1_sigma <- c(0, 0.1, 0.2)

head(primate.data.edited)
pd <- setNames(primate.data.edited$primate.tree2.tip.label, rownames(primate.data.edited))
pd

#### Next section ####

#Number of taxa to keep 
possible_taxa <- round(Ntip(primate.tree2)*c(.1, .5, 1))
str(possible_taxa)

covar_index <- length(possible_covar)
sigma_index <- length(test1_sigma)
replicate_index <- sequence(10)
covar_matrix <- matrix(test1_sigma[sigma_index], nrow=3, ncol=3)
covar_matrix[1,2] <- possible_covar[covar_index]
covar_matrix[2,1] <- covar_matrix[1,2]
traits <- sim.char(primate.tree2, nsim = 1, model = c("BM"), root=1)
sim.traits <- sim.char(primate.tree2, par = covar_matrix, 0.02, 100) ## Issues with this code, error in match.arg
sim.trait.test <- sim.char(primate.tree2, 0.2, 100) # Continuous character - unvariate 
sim.trait.fitC <- fitContinuous(primate.tree2, model = "BM", control = list(niter=10), ncores = 2)


taxa_index <- sequence(possible_taxa){
  pruned.tree <- geiger::drop.random(primate.tree2, n=ape::Ntip(primate.tree2) - possible_taxa[taxa_index])
  pruned.tree.all <- geiger::treedata(pruned.tree, traits)
  geiger.tree.results <- geiger::fitContinuous(pruned.tree.all$phy, dat = pruned.tree.all$data)
}

phy_pruned <- geiger::drop.random(primate.tree2, n=ape::Ntip(primate.tree2))
str(phy_pruned)
pruned_all <- geiger::treedata(phy_pruned, sim.trait.test)
geiger_results <- geiger::fitContinuous(phy_pruned$primate.tree2, dat=primate.data.edited)

local.results <- data.frame(
  covariation=possible_covar[covar_index],
  true.sigma=test1_sigma[sigma_index],
  ntip=possible_taxa[taxa_index],
  estimated.covar=
)
#### Next section ####

# 5. Figure out the likelihood 
BM.fit.lik <- fitContinuous(phy=primate.tree2, dat=primate.data.edited, SE=NA, control = list(niter=50), ncores = 2)
bm2 <- fitContinuous(phy=primate.tree2, dat=primate.data.edited, SE=NA, control = list(niter=50), ncores = 2)



#### Next section ####

# 6. Analyze the simulations 
#Create simulations

## Tests with out using the tree and data:
t<-0:100 #This sets the length of time for the simulation / number of generations 

sig2 <- 0.01 #This determines the relative size of jumps on the tree. When sigma squared is small, character changes tend to be smaller. Sigma squared has a maximum of 1, which suggests larger character changes. 

random.deviates.simulation <- rnorm(n=length(t) -1, sd=sqrt(sig2)) #This is the series of character state changes through time

stepwise_values <- numeric() # Generating trait values over time
stepwise_values[1] <-0 # Each trait value at time t is calculated by adding a change in trait value given in the vector to the trait value at t-1. Character state begins at 0. 
stepwise_values[2] <- stepwise_values[1] + random.deviates.simulation[1]
stepwise_values[3] <- stepwise_values[2] + random.deviates.simulation[2]
stepwise_values[4] <- stepwise_values[3] + random.deviates.simulation[3]
stepwise_values[5] <- stepwise_values[4] + random.deviates.simulation[4]


#Plot the trait values for the first 5 changes in time described above. 
abs_max <- max(abs(stepwise_values))
plot(1, stepwise_values[1], 
     xlim = c(0,5), ylim = c(-abs_max,abs_max),
     xlab = "t", ylab = "trait value",
     pch = 20, type = "o")

points(c(2:5), stepwise_values[2:5], pch=20)

lines(1:5, stepwise_values)

# Do this for all the values of x:
random.deviates.simulation <- c(0, cumsum(random.deviates.simulation))
plot(t, random.deviates.simulation, type = "l", ylim = c(-3, 3))

# Change the values of t and sigma squared:
t <- 0:100
sig2 <- 0.01
nsim <- 100
X <- matrix(rnorm(n = nsim * (length(t) - 1), sd = sqrt(sig2)), nsim, length(t) - 1)
sim_matrix <- cbind(rep(0, nsim), t(apply(X, 1, cumsum))) #matrix of simulations

# Plot simulations:
plot(t, sim_matrix[1, ], xlab = "time", ylab = "phenotype", ylim = c(-3, 3), type = "l")
apply(sim_matrix[2:nsim, ], 1, function(x, t) lines(t, x), t = t)
lines(t, sim_matrix[1, ], xlab = "time", ylab = "phenotype", ylim = c(-3, 3), col = "red") # the red is the first simulation run. 



#### Next section ####

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


##Projecting a continuous trait onto the branches of a tree using variable edge widths
hyoid.hornlengths <- setNames(data.primate1$data.T1, rownames(data.primate1))
object <- edge.widthMap(primate.tree2, hyoid.hornlengths)

#### Next section ####

# 8. Decide the minimum number of taxa needed/appropriate 


#### Next section ####

# 9. Compare results from two different primate tree 
# Comparison tree from 10k Trees Primate Project 
primate.tree3 <- ape::read.nexus("consensusTree_10kTrees_Primates_Version3_withCladeCredibilityValues.nex")
plot(primate.tree2)












































































