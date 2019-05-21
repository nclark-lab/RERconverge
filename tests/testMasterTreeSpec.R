#Try reading in yeast trees with and without specification of master tree
#Use first 19 trees to ensure there are not enough trees to build edge lengths
#Run from tests directory
source('../R/RERfuncs.R')
source('../R/RcppExports.R')
rightmaster = read.tree('yeast1.tre') #appropriate topology
wrongmaster = read.tree('yeast2.tre') #wrong topology
trtfail = readTrees('yeastfirst19.txt') #should fail (produce a master tree with no edge lengths)
trtsucceed = readTrees('yeastfirst19.txt',masterTree=rightmaster) #should use master tree provided
#Check error handling for incorrect topology of:
#a) master tree
trtfail2 = readTrees('yeastfirst19.txt',masterTree=wrongmaster)
#b) at least one gene tree
trtfail3 = readTrees('yeastfirst19_genetrees.txt',masterTree=rightmaster)
#Make sure RERs are calculable using the master tree
trer = getAllResiduals(trtsucceed)
#Error in apply(allbranch, 2, mean, na.rm = T, trim = mean.trim) :
#dim(X) must have a positive length
#Also try without rotating.... how does this work?
trtsucceed2 = readTrees('yeastfirst19.txt',masterTree=rightmaster,doRotate=F)
