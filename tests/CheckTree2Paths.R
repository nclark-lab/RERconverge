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
pdf('tests/herbivoreTest_al_wt_an_te_maxfgd.pdf')
par(mfrow=c(1,2))
herbv2state <- foreground2Tree(herb, toyTrees, clade="all", useSpecies=dietus)
herbv2weight <- foreground2Tree(herb, toyTrees, clade="weighted", useSpecies=dietus)
herbv2anc <- foreground2Tree(herb, toyTrees, clade="ancestral", useSpecies=dietus)
herbv2term <- foreground2Tree(herb, toyTrees, clade="terminal", useSpecies=dietus)
herbv2maxfgdclade <- foreground2Tree(herb, toyTrees, clade="MaximalForegroundClade", useSpecies=dietus)
dev.off()
