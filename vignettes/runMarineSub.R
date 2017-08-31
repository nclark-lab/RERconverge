if(F){
library("RERconverge")
#this will take some time
mamTrees=readTrees("../data/mammal_62_aa_sub.tre", max.read = 100)
#cutoff = exp(-7) based on the peak of the mean/variance plot


#this is the basic method
mamRER=getAllResiduals(mamTrees,useSpecies=mamTrees$masterTree$tip.label, transform = "none",weighted = F, cutoff=0.001)
#For some datasets scaling improves results
mamRERs=scale(mamRER)

#this method peforms better on benchmarks
mamRERlogW=getAllResiduals(mamTrees,useSpecies=mamTrees$masterTree$tip.label, transform = "log",weighted = T, cutoff=0.001)

#again can be scaled, almost always improves results in our exprience
mamRERlogWs=scale(mamRERlogW)

#read the binary tree
marineb=read.tree("../data/MarineTreeBin.txt")
plot(marineb)

#we can also do this from the foreground set but only the extant branches will be set to 1
foreground=c("triMan1", "turTru2", "odoRosiDi", "lepWed1")
phenvMarine2=foreground2Paths(foreground, mamTrees)


#convert it to a paths vector
  phenvMarine=tree2Paths(marineb, mamTrees)



corMarine=getAllCor(mamRER, phenvMarine)
hist(corMarine$P)

corMarineLogW=getAllCor(mamRERlogW, phenvMarine)
hist(corMarineLogW$P)

corMarineLogWs=getAllCor(mamRERlogWs, phenvMarine)
hist(corMarineLogWs$P)
}
