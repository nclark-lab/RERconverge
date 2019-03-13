#For testing new functions to plot and save RERs by gene
#(new plotting functions in plottingFuncs.R and tree saving functions in RERfuncs.R)
if (!require("RERconverge", character.only=T, quietly=T)) {
  require(devtools)
  install_github("nclark-lab/RERconverge",
                 ref="AddressReviewerComments") #can be modified to specify a particular branch
}
library(RERconverge)
#Source the functions locally for testing
repodir='~/repos/RERconverge' #replace with local directory
source(paste(repodir,'/R/plottingFuncs.R',sep=''))
data("toyTrees")
data("mamRERw")
phenvExample <- foreground2Paths(c("Vole","Squirrel"),toyTrees,clade="terminal")
relGene = "BEND3"
#Find a way to map RERs to tree edges in order to use treePlotNew