if(T){
  if (!require("RERconverge", character.only=T, quietly=T)) {
    require(devtools)
    install_github("nclark-lab/RERconverge")
  }
#the following will take some time
mamTrees=readTrees("../data/mammal62aa_meredplus_wCM.trees", max.read = 200)
#cutoff = exp(-7) based on the peak of the mean/variance plot


#this is the basic method (not demonstrated in the vignette)
mamRER=getAllResiduals(mamTrees,useSpecies=mamTrees$masterTree$tip.label, transform = "none",weighted = F, scale = T)
#For some datasets scaling improves results

#this method peforms better on benchmarks, weighting works with all transforms but sqrt is recommended
#(scaling no longer recommended and not demonstrated in the vignette)
mamRERw=getAllResiduals(mamTrees,useSpecies=mamTrees$masterTree$tip.label, transform = "sqrt",weighted = T, scale=T)

#this method peforms better on benchmarks, weighting works with all transforms but sqrt is recommended
mamRERv=getAllResiduals(mamTrees,useSpecies=mamTrees$masterTree$tip.label,
                       transform = "sqrt", weighted = T, cutoff=0.001)




#read the binary tree
marineb=read.tree("data/MarineTreeBin.txt")
#offset 0 length edges, makes no difference for further computation
marineb$edge.length[marineb$edge.length<1]=0.05
plot(marineb)

#we can also do this from the foreground set but only the extant branches will be set to 1
foreground=c("triMan1", "orcOrc1", "turTru2", "odoRosDi", "lepWed1")
plot(marineb)

marineb2=foreground2Tree(foreground, mamTrees)
marineb3=foreground2Tree(foreground, mamTrees, collapse2anc = F)


#convert a marine tree to a paths vector
phenvMarine=tree2Paths(marineb, mamTrees)
phenvMarine

corMarine = getAllCor(mamRERv, phenvMarine)


#corMarine=getAllCor(mamRER, phenvMarine)
#hist(corMarine$P)

#corMarineW=getAllCor(mamRERw, phenvMarine)
#hist(corMarineW$P)
#qqplot(corMarine$P, corMarineW$P, log="xy");abline(a=0,b=1)

}
