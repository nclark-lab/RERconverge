#Load toy data and use to test method for estimating difference in RERs
require(phytools)
require(RERconverge)
#source RERfuncs locally
source('../R/RERfuncs.R')
#need to re-do readTrees to get reoriented trees
noneutherians <- c("Platypus","Wallaby","Tasmanian_devil","Opossum")
toyTreesNew = readTrees('../ext/subsetMammalGeneTrees.txt',reorient=T,outgroup=noneutherians)
#latest pseudo-rooting problem is for a tree missing platypus
#how will this tree/path ultimately be stored in the trees object?

#convert to RER trees
rertrees = returnRersAsTreesAll(toyTrees,mamRERw)
#multiple upstream branches
diffrertrees = list()
for (i in 1:length(rertrees)) {
  print(i)
  diffrertrees[[i]] = diffBranches(rertrees[[i]])
}


ttree = toyTrees$masterTree
#re-root tree
rttree = midpoint.root(ttree)
#plot these to visualize
par(mfrow=c(1,2))
plot(ttree,cex=0.5,no.margin=T)
nodelabels(cex=0.5)
plot(rttree,cex=0.5,no.margin=T)
nodelabels(cex=0.5)

#map new edges to original edges
#edge X should be edge Y - edge Z
#worst case: human
wtt = which(ttree$tip.label=="Human") #62
wte = which(ttree$edge[,2] == wtt) #121
#want to replace this with (internal branch - terminal branch) from rooted tree

#make example tree
egspecies = c("Human","Chimp","Gorilla","Cow","Platypus")
ttsub = drop.tip(ttree,ttree$tip.label[ttree$tip.label %in% egspecies == F])
rttsub = drop.tip(rttree,rttree$tip.label[rttree$tip.label %in% egspecies == F])
par(mfrow=c(1,2))
plot(ttsub,cex=0.5,no.margin=T)
nodelabels(cex=0.5)
plot(rttsub,cex=0.5,no.margin=T)
nodelabels(cex=0.5)

#switch to test distances
newrttsub=read.tree(text="((((Chimp:1,Human:2):100,Gorilla:200):1000,Cow:2000):10000,Platypus:20000);")
newttsub = root(unroot(newrttsub),"Human")
write.tree(newttsub,file="../tests/testTree.tre")

#Maria has a full tree and one that is missing the non-placentals
#Goal: revise readTrees so that the master tree is rooted at midpoint root
#topology of all sub-trees should match that of the master tree
#then the edges should be in the "right direction" for subtraction
