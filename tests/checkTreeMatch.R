#Check that RERconverge will identify issues with tree topology,
#namely, ones where the trait tree is not a subset of the master gene tree.
require(devtools)
#install_github("nclark-lab/RERconverge",branch="MatchingTrees")
library('RERconverge')
#repopath = '~/repos/RERconverge' #local path to repo
repopath = '~/Documents/GitHub/RERconverge'
rerpath = find.package('RERconverge') #will files in the tests directory be installed here?

#Use same master tree for all cases:
toytreefile = "subsetMammalGeneTrees.txt"
toyTrees=readTrees(paste(rerpath,"/extdata/",toytreefile,sep=""), max.read = 200)

#Example cases:
# -rooted trait tree
# -different tree topology
#need a more recent ucsc tree.... this one has numbers in it rather than species names
ucsctree = read.tree(paste(repopath,'tests/ucsc62mammals_commonnames.tre',sep='/'))
ttucsc = tree2Paths(ucsctree, toyTrees)

# -species in the trait tree that are not in the master tree
pon1tree = read.tree(paste(repopath,'tests/PON1_addl_species.tre',sep='/'))
ttpon1 = tree2Paths(pon1tree, toyTrees)

# -different species filtered from trait tree and master tree
