#Load toy data and use to test method for estimating difference in RERs
require(phytools)
require(RERconverge)
#source RERfuncs locally
source('../R/RERfuncs.R')
source('../R/RcppExports.R')
#need to re-do readTrees to get reoriented trees
noneutherians <- c("Platypus","Wallaby","Tasmanian_devil","Opossum")
toyTreesNew = readTrees('../ext/subsetMammalGeneTrees.txt',reorient=T,outgroup=noneutherians)
#It works?

#estimate RERs
ttResid = getAllResiduals(toyTreesNew)

#The below two steps are now combined in changeInRers in RERfuncs.
#convert to RER trees
rertrees = returnRersAsTreesAll(toyTreesNew,ttResid)

#get delta_RER trees and convert back to RER matrix
diffrertrees = list()
diffrermat = matrix(data=NA,nrow=nrow(ttResid),ncol=ncol(ttResid))
colnames(diffrermat) = colnames(ttResid)
rownames(diffrermat) = rownames(ttResid)
for (i in 1:length(rertrees)) {
  #print(i)
  diffrertrees[[i]] = diffBranches(rertrees[[i]])
  ee=edgeIndexRelativeMaster(diffrertrees[[i]], toyTreesNew$masterTree)
  #this is the columns of the  paths matrix
  ii= toyTreesNew$matIndex[ee[, c(2,1)]]
  diffrermat[i,ii] <- diffrertrees[[i]]$edge.length
}

#test whether this can be plotted
phenvExample <- foreground2Paths(c("Vole","Squirrel"),toyTreesNew,clade="terminal")
plotRers(diffrermat,"BEND3",phenv=phenvExample)


