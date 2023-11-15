#This demonstrates the issues with inferring a binary tree
#using clade="all"

#Run from RERconverge repo main directory as working directory

#Source RERfuncs from MatchingTrees branch:
library(RERconverge)

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
pdf('tests/herbivoreTest_allOptions.pdf')
par(mfrow=c(2,2))
herbv2statebi <- foreground2Tree(herb, toyTrees, clade="all", transition = 'bidirectional',useSpecies=dietus)
herbv2stateuni <- foreground2Tree(herb, toyTrees, clade="all", transition = 'unidirectional', useSpecies=dietus)
herbv2statebiwt <- foreground2Tree(herb, toyTrees, clade="all", transition = 'bidirectional', weighted = T, useSpecies=dietus)
herbv2stateuniwt <- foreground2Tree(herb, toyTrees, clade="all", transition = 'unidirectional', weighted = T, useSpecies=dietus)
dev.off()
pdf('tests/herbivoreTest_anctermOptions.pdf')
par(mfrow=c(2,2))
herbv2ancbi <- foreground2Tree(herb, toyTrees, clade="ancestral", transition = 'bidirectional',useSpecies=dietus)
herbv2ancuni <- foreground2Tree(herb, toyTrees, clade="ancestral", transition = 'unidirectional', useSpecies=dietus)
herbv2term <- foreground2Tree(herb, toyTrees, clade="terminal", useSpecies=dietus)
herbv2termwt <- foreground2Tree(herb, toyTrees, clade="terminal", weighted = T, useSpecies=dietus)
dev.off()
