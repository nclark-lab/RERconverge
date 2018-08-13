#Early version of the steps used in the walk-through
#If including, modify script to align with latest vignette

library("RERconverge")
rerpath = paste(.libPaths()[1],"/RERconverge",sep="")
toytreefile = "mammal62aa_meredplus_wCM.trees" #change filename once toy trees available
#this will take some time
mamTrees=readTrees(paste(rerpath,"/extdata/",toytreefile,sep=""), max.read = 200)


#this is the basic method (not used in the vignette)
mamRER=getAllResiduals(mamTrees,useSpecies=mamTrees$masterTree$tip.label, 
                       transform = "none",weighted = F, scale = T)


#this method peforms better on benchmarks and is used in the vignette
#weighting works with all transforms but sqrt is recommended
mamRERw=getAllResiduals(mamTrees,useSpecies=mamTrees$masterTree$tip.label, 
                        transform = "sqrt",weighted = T, scale=T)


#read the binary tree
marineb=read.tree(paste(rerpath,"/extdata/MarineTreeBin.txt",sep=""))
#offset 0 length edges, makes no difference for further computation (double check this)
marineb$edge.length[marineb$edge.length<1]=0.05
plot(marineb)

#we can also do this from the foreground set
foreground=c("triMan1", "orcOrc1", "turTru2", "odoRosDi", "lepWed1")
#there are two options, the default will put all the weight for a clade on the ancestral branch
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



