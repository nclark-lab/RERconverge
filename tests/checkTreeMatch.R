#Check that RERconverge will identify issues with tree topology,
#namely, ones where the trait tree is not a subset of the master gene tree.
library('RERconverge')
rerpath = find.package('RERconverge')

#Use same master tree for all cases:
toytreefile = "subsetMammalGeneTrees.txt" 
toyTrees=readTrees(paste(rerpath,"/extdata/",toytreefile,sep=""), max.read = 200)

#Example cases:
# -rooted trait tree
# -different tree topology
# -species in the trait tree that are not in the master tree
pon1tree = read.tree(paste(rerpath,'PON1_addlspecies.tre',sep='/'))
ttpon1 = tree2Paths(pon1tree, toyTrees)

# -different species filtered from trait tree and master tree
