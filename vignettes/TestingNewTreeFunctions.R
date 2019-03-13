#For testing new functions to plot and save RERs by gene
#(new plotting functions in plottingFuncs.R and tree saving functions in RERfuncs.R)
if (!require("RERconverge", character.only=T, quietly=T)) {
  require(devtools)
  install_github("nclark-lab/RERconverge",
                 ref="AddressReviewerComments") #can be modified to specify a particular branch
}
library(RERconverge)
rerpath = find.package('RERconverge')
toytreefile = "subsetMammalGeneTrees.txt" 
toyTrees=readTrees(paste(rerpath,"/extdata/",toytreefile,sep=""), max.read = 200)
data("logAdultWeightcm") 
mamRERw = getAllResiduals(toyTrees,useSpecies=names(logAdultWeightcm), 
                          transform = "sqrt", weighted = T, scale = T)
phenvExample <- foreground2Paths(c("Vole","Squirrel"),toyTrees,clade="terminal")
