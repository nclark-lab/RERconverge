#add to readTrees after rotateConstr to "reorient" the trees

#root the master tree
mytreeFull.rooted=midpoint.root(mytreeFull)

#root a new tree with missing species
#
#prunt the rooted masterTree to the correct species
prunedTree=pruneTree(mytreeFull.rooted, mytreeMissing$tip.label)
#find the new root node
rootNode=setdiff(prunedTree$edge[,1], prunedTree$edge[,2])
#find the first child or the root node
rootNode1down=prunedTree$edge[prunedTree$edge[,1]==rootNode,2][1]
#get the corresponding outgroup
outgroup=treeTraverse(prunedTree, rootNode1down)
#root the tree with missing species on the outgroup
mytreeMissing.rooted=root(mytreeMissing, outgroup = outgroup)

plot(mytreeMissing.rooted)

#at the end these need to be unrooted again