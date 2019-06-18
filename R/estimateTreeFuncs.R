#add-on functions for creating trees for use with RERconverge
require(phangorn)
require(phytools)
require(tools)
pruneAlnFromTree = function(alnfile, treefile, type = "AA", format = "fasta", writealn=TRUE) {
  #prune the alignment to have only the species in the tree
  #read in the alignment
  alnPhyDat = read.phyDat(alnfile, type = type, format = format)
  #read in the treefile
  genetree = read.tree(treefile)
  #eliminate species in the alignment but not the tree
  inboth = intersect(names(alnPhyDat),genetree$tip.label)
  alnPhyDat = subset(alnPhyDat, subset = inboth)
  if (writealn) {
    #write the new alignment with a revised filename
    fe = file_ext(alnfile)
    fpse = file_path_sans_ext(alnfile)
    write.phyDat(alnPhyDat, file=paste(fpse,".pruned.",fe,sep=""), format = format)
  }
  return(alnPhyDat)
}
pruneTreeFromAln = function (treefile, alnfile, type = "AA", format = "fasta", writetree=TRUE) {
  #prune the tree to have only the species in the alignment
  #read in the alignment
  alnPhyDat = read.phyDat(alnfile, type = type, format = format)
  #read in the treefile
  genetree = read.tree(treefile)
  #eliminate species in the alignment but not the tree and vice versa
  inboth = intersect(names(alnPhyDat),genetree$tip.label)
  todrop = genetree$tip.label[genetree$tip.label %in% inboth == FALSE]
  if (length(todrop) > 0) {
    genetree = drop.tip(genetree, todrop)
  }
  #unroot the tree
  genetree = unroot(genetree)
  if (writetree) {
    #write the new tree with a revised filename
    fe = file_ext(treefile)
    fpse = file_path_sans_ext(treefile)
    write.tree(genetree, file=paste(fpse,".pruned.",fe,sep=""))
  }
  return(genetree)
}

#' Estimate a ML tree from a given topology
#' Uses `pml` and `optim.pml` from the `phangorn` package to estimate the tree.

#' @param alnfile The path to the sequence alignment file
#' @param treefile The path to the master tree file (whose topology will be used to generate the tree)
#' @param submodel Substitution model to be used to estimate the tree
#' @param type "AA" for amino acid or "DNA" for DNA
#' @param format Format of the alignment file (commonly used formats include fasta and phylip)
#' @param k Number of intervals in the discrete gamma distribution for `pml`
#' @param ... Further arguments passed to `pml` or `optim.pml`
#' @return A list: `tree.opt` is the optimized output of `optim.pml`; 
#' `tree.init` is the initial tree estimated by `pml`
#' @seealso \code{\link[phangorn]{phyDat}} for alignment formats,
#'  \code{\link[phangorn]{pml}} and \code{\link[phangorn]{optim.pml}} for tree estimation
#' @export
estimatePhangornTree = function(alnfile, treefile, submodel="LG", type = "AA", 
                                format = "fasta", k=4, ...) {
  #Generates distance-based tree using submodel for substitutions
  #Read in pruned alignment and pruned master tree
  #If species lists are different, prune both to be the same.
  #read in the alignment
  alnPhyDat = read.phyDat(alnfile, type = type, format = format)
  #read in the treefile
  genetree = read.tree(treefile)
  #eliminate species in the alignment but not the tree and vice versa
  inboth = intersect(names(alnPhyDat),genetree$tip.label)
  todropg = genetree$tip.label[genetree$tip.label %in% inboth == FALSE]
  if (length(todropg) > 0) {
    genetree = drop.tip(genetree, todropg)
  }
  if (length(inboth) < length(names(alnPhyDat))) {
    alnPhyDat = subset(alnPhyDat, subset = inboth)
  }
  #unroot the tree
  genetree = unroot(genetree)
  #Run distance estimation using submodel
  #generate an initial pml tree
  lgptree = pml(genetree, alnPhyDat, model = submodel, k = k, rearrangement="none", ...)
  #generate a tree
  lgopttree = optim.pml(lgptree,optInv=T,optGamma=T,optEdge=T,rearrangement="none",
                        model=submodel, ...)
  return(list("tree.init"=lgptree,"tree.opt"=lgopttree))
}