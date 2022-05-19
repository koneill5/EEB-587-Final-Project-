# EEB 587 Final Project - due the last day of finals (18 May 2022)

# 1. Create project question 
#Power analysis to decide the appropriate number/range of taxa needed for dissertation work 

# Question 1 - What is an appropriate number of taxa to use when estimating the co-variation of two traits?

# Question 2 - Does this minimum number of taxa differ when using two different primate phylogenies? 
 # source for the whole thing so it runs and fix some of the non-need items
#Methods/steps:

#2. Input and read primate tree 
library(ape)
library(geiger)
library(phytools)
library(Runuran)
library(tidyverse)
library(ggplot2)
library(ggtree)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggtree")
library(dplyr)

#Use the following tree:
#Downloaded from the 10k Trees Primates program 
#Branch lengths are proportional to absolute time. Ultrametric tree as Chronogram 
primate.tree2 <- ape::read.nexus("consensusTree_10kTrees_Primates_Version3.nex")
plot(primate.tree2, show.tip.label =FALSE)
primate.tree2$edge
str(primate.tree2) #structure of the tree
primate.tree2$tip.label #vector of names 

is.binary(primate.tree2)
is.ultrametric(primate.tree2)
is.rooted(primate.tree2)
multi2di(primate.tree2)

better_tree <- ape::multi2di(primate.tree2)
is.binary(better_tree)

#### Next Section ####

# 3. Simulate data
# Want to estimate co-variation between two traits 
# Possible co-variation values:
possible_covar <- c(0, 0.5, 1)

#Trait 1 (Long-axis length of thyrohyal [greater horn] of hyoid) average: 10 mm
Trait_1 <- 10

#Trait 2 (Length of clavicle) average: 145 mm 
Trait_2 <- 145

# sim trait on geiger is a better way to simulate data 

# Generate random data:
set.seed(1234)
x <- rnorm(n=301, mean=10, sd=5)
require("Runuran")
data.T1 <- urnorm(n=301, mean=10, sd=2,lb=1, ub=100)
data.T2 <- urnorm(n=301, mean= 145, sd=25,lb=1, ub=200)
dt <- cbind(data.T1, data.T2)
data.primate <- cbind(primate.tree2$tip.label, data.T1, data.T2) # set as rownames()
data.primate1 <- as.data.frame(data.primate)
rownames(data.primate1) = data.primate1$V1
data.primate1
head(data.primate1)


# Compare taxa in dat and tree:
treedata(primate.tree2, data.primate1)
name.check(primate.tree2,data.primate1)

#### Next section ####

# 4. Figure out a good sigma^2 
test1_sigma <- c(0.01, 0.05)


#### Next section ####

#Number of taxa to keep 
possible_taxa <- round(Ntip(primate.tree2)*c(.1, .5, 1))
str(possible_taxa)

covar_index <- length(possible_covar)
sigma_index <- length(test1_sigma)
replicate_index <- sequence(10)
covar_matrix <- matrix(test1_sigma[sigma_index], nrow=2, ncol=2)
covar_matrix[1,2] <- possible_covar[covar_index]
covar_matrix[2,1] <- covar_matrix[1,2]
traits <- sim.char(primate.tree2, nsim = 1, model = c("BM"), root=1)
sim.traits <- sim.char(primate.tree2, par = covar_matrix, model= "BM", nsim= 100)## Issues with this code, error in match.arg / starting at 0 by default 
head(sim.traits)
plot(sim.traits)

sim.traits <- setNames(sim.traits, rownames(sim.traits))

sim.2 <- sim.char(better_tree, par = covar_matrix, model = "BM", nsim = 301)
plot(sim.2)



# fitContinuous tries:
head(sort(names(better_tree)))
sim.trait.fitC <- fitContinuous(phy=better_tree,dat=data.primate1$data.T1,SE=NA, control = list(niter=10), ncores = 2)# error 'phy' is not a binary tree
class(data.primate1$data.T1)
data.primate1$data.T1 = as.integer(data.primate1$data.T1)
data.primate1$data.T2 = as.integer(data.primate1$data.T2)
class(data.primate1$data.T1)
class(data.primate1$data.T2)
sim.fitC <- fitContinuous(ape::multi2di(primate.tree2),data.primate1$data.T1)

all.things = treedata(better_tree,data.primate1)
phy.2 <- all.things$phy
dat.2 <- all.things$data
fittest <- fitContinuous(phy.2, dat.2[,"data.T1"], SE=NA, control = list(niter=50), ncores = 2)


  ## Work on this section: 
taxa_index <- seq(possible_taxa){
  pruned.tree <- geiger::drop.random(primate.tree2, n=ape::Ntip(primate.tree2) - possible_taxa[taxa_index])
  pruned.tree.all <- geiger::treedata(pruned.tree, traits)
  geiger.tree.results <- geiger::fitContinuous(pruned.tree.all$phy, dat = pruned.tree.all$data)
}


phy_pruned <- geiger::drop.random(primate.tree2, n=ape::Ntip(primate.tree2) - possible_taxa)
str(phy_pruned)
phy_pruned$tip.label
treedata(better_tree, sim.traits)

pruned_all <- geiger::treedata(phy_pruned, sim.traits,warnings = TRUE) #error for incorrect dimensions, need to subset the specific taxa named in pruned? 
geiger_results_test <- geiger::fitContinuous(phy_pruned$primate.tree2, dat=data.primate1)

results <- data.frame()
for (covar_index in seq(length(possible_covar))){
  for (sigma_index in seq(length(test1_sigma))) {
    for (replicate_index in seq(10)) {
      covar_matrix <- matrix(test1_sigma[sigma_index], nrow=2, ncol=2)
      covar_matrix[1,2] <- possible_covar[covar_index]
      covar_matrix[2,1] <- covar_matrix[1,2]
      sim.traits <- sim.char(better_tree, par = covar_matrix, model= "BM", nsim= 100)
      for (taxa_index in seq(possible_taxa)) {
        geiger_results<- geiger::fitContinuous(better_tree,sim.traits)
      }
    }
  }
}

local.results <- data.frame(
  covariation=possible_covar[covar_index],
  true.sigma=test1_sigma[sigma_index],
  ntip=possible_taxa[taxa_index],
  estimated.covar=geiger_results_test$opt
  estimated_sigma_squared=geiger_results_test$opt
  replicate_number=replicate_index
)

#### Next section ####

# 5. Figure out the likelihood and # 6. Analyze the simulations 
BM.fit.lik <- fitContinuous(phy=better_tree, dat=hl)
bm2 <- fitContinuous(phy=better_tree, dat=hl, SE=NA, control = list(niter=50), ncores = 2)


#Computer phylogenetic signal with two methods:

all.things = treedata(better_tree,data.primate1)
phy.2 <- all.things$phy
dat.2 <- all.things$data

hl <- setNames(data.primate1$data.T1, rownames(data.primate1))
cl <- setNames(data.primate1$data.T2, rownames(data.primate1))

k.hl <- phylosig(better_tree,hl,method = "K", test = FALSE, nsim=1000)
print(k.hl)
plot(k.hl)

k.cl <- phylosig(better_tree,cl, method = "K", test = FALSE, nsim=1000)
print(k.cl)
plot(k.cl)

l.hl <- phylosig(better_tree, hl, method = "lambda", test = TRUE)
print(l.hl)
plot(l.hl)

l.cl <- phylosig(better_tree,cl,method = "lambda", test = TRUE)
print(l.cl)
plot(l.cl)

fct.hl <- fitContinuous(better_tree,hl)
print(names(fct.hl))
print(fct.hl)
fct.hl$opt


fct.cl <- fitContinuous(better_tree,cl)
print(names(fct.cl))
print(fct.cl)
fct.cl$opt

flik.hl=fct.hl$lik
print(argn(flik.hl))

flik.cl=fct.cl$lik
print(argn(flik.cl))


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


# 7. Visualize 
## box plots are okay 
N<-50
t <- fastBM(better_tree,nsim = 10)
fn <- function(better_tree,x,sampling.var){
  foo<-function()rep(sampling.var,N)/sample(1:5,size=N,replace=TRUE)
  v<-replicate(nsim,foo(),simplify=FALSE)
  foo<-function(x,v) mapply(sampleFrom,xbar=x,xvar=v)
  mapply(foo,x=x,v=v,SIMPLIFY=FALSE)
}

ve <- seq(0,0.4,by=0.05) #this is going to test different sigma2
K <- lambda <- matrix(NA,nsim,length(ve),dimnames = list(NULL,ve))
xe <- list()
for (i in 1:length(ve)) {
  xe[[i]]<-fn(better_tree,hl,ve[i])
  K[,i]<-mapply(phylosig,tree=better_tree,x=xe[[i]])
  foo<-function(better_tree,hl) phylosig(better_tree,hl,method = "lambda")$lambda
  lambda[,i]<-mapply(foo,tree=better_tree,x=xe[[i]])
}

colMeans(K)
boxplot(K, xlab="Sample size (number of taxa)" , ylab="sigma2 variance")

boxplot(horn.length)

##Projecting a continuous trait onto the branches of a tree using variable edge widths
#Map continuous traits on the tree 
horn.length <- setNames(data.primate1$data.T1, row.names(data.primate1))
object.plot <- edge.widthMap(better_tree,horn.length, min.width=0.05)
plot(object.plot, border=1.5, color=palette()[2], legend = "horn length")

clav.length <- setNames(data.primate1$data.T2, rownames(data.primate1))
obj.plot2 <- edge.widthMap(better_tree,clav.length, min.width=0.05)
plot(obj.plot2, border=1.5, color=palette()[5], legend = "clavicle length")

#### Next section ####
#Contmap 
primate.contmap.horn <- contMap(better_tree, horn.length)
plot(primate.contmap.horn, fsize=c(0.7,0.8), leg.txt="horn length")
par(mar=c(5.1,4.1,4.1,2.1))
primate.contmap.clav <- contMap(better_tree,clav.length)
plot(primate.contmap.clav, fsize=c(0.7,0.8), leg.txt="clavicle length")
par(mar=c(5.1,4.1,4.1,2.1))

all.things = treedata(better_tree,data.primate1)
phy.2 <- all.things$phy
dat.2 <- all.things$data

test.frame <- setNames(all.things$phy,all.things$data)

par(mfrow=c(1,2))
plotTree(better_tree,mar=c(5.1,1.1,2.1,0.1))
par(mar=c(5.1,0.1,2.1,1.1))
boxplot(horn.length~factor(names(horn.length),levels = better_tree$tip.label), horizontal = TRUE, axes=FALSE, xlim=c(1,Ntip(better_tree)))
axis(1)
title(xlab = "horn length")

par(mfrow=c(1,2))
plotTree(better_tree,mar=c(5.1,1.1,2.1,0.1))
par(mar=c(5.1,0.1,2.1,1.1))
boxplot(clav.length~factor(names(clav.length),levels = better_tree$tip.label), horizontal = TRUE, axes=FALSE, xlim=c(1,Ntip(better_tree)))
axis(1)
title(xlab = "clavicle length")



data.primate1
op <- par(no.readonly = TRUE)
par(mfrow=c(1,2))
with(data.primate1, plot(data.primate1$V1,data.primate1$data.T1))

# 8. Decide the minimum number of taxa needed/appropriate 


#### Next section ####

# 9. Compare results from two different primate tree 
# Comparison tree from 10k Trees Primate Project 
primate.tree3 <- ape::read.nexus("consensusTree_10kTrees_Primates_Version3_withCladeCredibilityValues.nex")
plot(primate.tree2)
primate.tree3$tip.label
primate.tree3$edge
str(primate.tree3)









































































