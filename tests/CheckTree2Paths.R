#This demonstrates the issues with inferring a binary tree
#using clade="all"

#Run from RERconverge repo main directory as working directory

#Source RERfuncs from MatchingTrees branch:
source('R/RERfuncs.R')

#Load toy trees data:
rerpath = find.package('RERconverge')
toytreefile = "subsetMammalGeneTrees.txt"
toyTrees=readTrees(paste(rerpath,"/extdata/",toytreefile,sep=""), max.read = 200)

#Herbivore data:
load('tests/herbivoreTest.RData')

#Dietary categories
print('Herbivore species:')
print(herb)
print('Non-herbivore species:')
print(dietus[dietus %in% herb == F])

#Plot output:
pdf('tests/herbivoreTest.pdf')
herbv2state <- foreground2Tree(herb, toyTrees, clade="all", useSpecies=dietus)
dev.off()
