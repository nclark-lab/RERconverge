
library("RERconverge")
#this will take some time
mamTrees=readTrees("data/mammal_62_aa_sub.tre", max.read = 1000)
#cutoff = exp(-7) based on the peak of the mean/variance plot


#this is the basic method
mamRER=getAllResiduals(mamTrees,useSpecies=mamTrees$masterTree$tip.label, transform = "none",weighted = F, scale = T)


#this method peforms better on benchmarks, weighting works with all transforms but sqrt is recommended
mamRERw=getAllResiduals(mamTrees,useSpecies=mamTrees$masterTree$tip.label, transform = "log",weighted = T, scale=T)




#read the binary tree
marineb=read.tree("data/MarineTreeBin.txt")
#offset 0 length edges, makes no difference for further computation
marineb$edge.length[marineb$edge.length<1]=0.05
plot(marineb)

#we can also do this from the foreground set
foreground=c("triMan1", "orcOrc1", "turTru2", "odoRosDi", "lepWed1")
#there are two options, the default will put all the weight for a clade on the ancestral branch
par(mfrow=c(1,3))
marineb1=foreground2Tree(foreground, mamTrees, collapse2anc = T)
#alternatively the weight will be put on the extant branch but this approach suffers from phylogenetic dependence
marineb2=foreground2Tree(foreground, mamTrees, collapse2anc = F)
#another approach is to spread the weight along the entire clade, like this. In order to handle the phylogenetic dependence in this case we case use weighted correlation, see bellow
marineb3=foreground2Tree(foreground, mamTrees, wholeClade = T)



#convert a marine tree to a paths vector
phenvMarine=tree2Paths(marineb1, mamTrees)
#do the correlation
corMarine1=getAllCor(mamRERw, phenvMarine)

#we can also do this with the other trees we generated
corMarine2=getAllCor(mamRERw, tree2Paths(marineb2, mamTrees))

#If we want to consider the whole clade as a single observation (using tree version 3) we need to use weighted correlation
corMarine3=getAllCor(mamRERw, tree2Paths(marineb3, mamTrees), weights = T)



