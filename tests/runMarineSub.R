




library("RERconverge")
#this will take some time
marineTrees=readTrees("ext/mammal_62_aa_sub.tre", max.read = 1000)



#run geAllREsiduals with default parameters
RERw=getAllResiduals(marineTrees)



#If we have a foreground set of species we can build a tree from it

foreground=c("triMan1", "orcOrc1", "turTru2", "odoRosDi", "lepWed1")
#but if the foreground species come in clades there are mutliple ways to do this, we will explore all of them
par(mfrow=c(2,2))
#naively we can set *all* the branch leading up to a forgreound species to 1
marinebAll=foreground2Tree(foreground, marineTrees, plot=T, clade="all");title("all")
#alternatively we can use just the ancestral branches
marinebAnc=foreground2Tree(foreground, marineTrees, clade="ancestral");title("ancestral")
#we can also use just the terminal branches
marinebTerm=foreground2Tree(foreground, marineTrees, clade="terminal");title("terminal")
#Finally, we can treat each clade as a single observation which lets us use all the branches but downweight them for the  phylogenetic dependence. This case is also handles differently when correlations are computed
marinebW=foreground2Tree(foreground, marineTrees, clade="weighted")
#note that the triTur and orcOrc clade is now weighted 1/3

corMarineAll=correlateWithBinaryPhenotype(RERw, tree2Paths(marinebAll, marineTrees))
corMarineAnc=correlateWithBinaryPhenotype(RERw, tree2Paths(marinebAnc, marineTrees))
corMarineTerm=correlateWithBinaryPhenotype(RERw, tree2Paths(marinebTerm, marineTrees))
corMarineW=correlateWithBinaryPhenotype(RERw, tree2Paths(marinebW, marineTrees))













