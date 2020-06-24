#add-on functions for creating trees for use with RERconverge
require(phangorn)
require(phytools)
require(tools) #may be unnecessary for tree estimation functions but useful for
#versioning in R
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
#' @return A list:
#' `tree.opt` is the tree from the optimized output of `optim.pml`;
#' `results.opt` is the optimized output of `optim.pml`;
#' `results.init` is the initial results estimated by `pml`
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
  #just in case, set all branches to 1 first (pml abhors a vacuum... or a zero)
  genetree$edge.length = c(rep(1,length(genetree$edge.length)))
  #Run distance estimation using submodel
  #generate an initial pml tree
  lgptree = pml(genetree, alnPhyDat, model = submodel, k = k, rearrangement="none", ...)
  #generate a tree
  #use capture.output to suppress optimization output?
  lgopttree = optim.pml(lgptree,optInv=T,optGamma=T,optEdge=T,rearrangement="none",
                        model=submodel, ...)
  return(list("results.init"=lgptree,"results.opt"=lgopttree, "tree.opt"=lgopttree$tree))
}


#' Estimates a ML trees from a given topology for a list of alignments. Wirtes a text file compatble with readTrees().
#' Uses `pml` and `optim.pml` from the `phangorn` package to estimate the tree.

#' @param alndir The path to directories with alignment files. Trees will be named with the alignment file name after stripping the extension. Either alnfiles or alndir must be supplied.
#' @param alnfiles A character vector of paths to alignment files. Such as one produced by \code{\link{list.files}}.
#' Either alnfiles or alndir must be supplied.
#' @param pattern An optional regular expression for files in the alndir director. As in  ".*fasta".
#' @param treefile The path to the master tree file (whose topology will be used to generate the tree)
#' @param output.file The file where the output trees will be written. This file can be read with readTrees().
#' @param submodel Substitution model to be used to estimate the tree
#' @param type "AA" for amino acid or "DNA" for DNA
#' @param format Format of the alignment file (commonly used formats include fasta and phylip)
#' @param k Number of intervals in the discrete gamma distribution for `pml`
#' @param ... Further arguments passed to `pml` or `optim.pml`
#' @return A list:
#' `tree.opt` is the tree from the optimized output of `optim.pml`;
#' `results.opt` is the optimized output of `optim.pml`;
#' `results.init` is the initial results estimated by `pml`
#' @seealso \code{\link[phangorn]{phyDat}} for alignment formats,
#'  \code{\link[phangorn]{pml}} and \code{\link[phangorn]{optim.pml}} for tree estimation
#' @export
estimatePhangornTreeAll = function(alnfiles=NULL, alndir=NULL, pattern=NULL, treefile, output.file=NULL,
                                   submodel="LG", type = "AA", format = "fasta", k=4, ...) {
  if(is.null(output.file)){
    stop("output.file must be supplied")
  }
  if (is.null(alnfiles)&&is.null(alndir)){
    stop("Either alnfiles or alndir must be supplied")
  }

  if (!is.null(alnfiles)&&!is.null(alndir)){
    stop("Only one of alnfiles or alndir must be supplied")
  }

  if (!is.null(alndir)){
    alnfiles=list.files(path=alndir, pattern=pattern, full.names = TRUE)
  }
  names=basename(alnfiles)
  names=sub("\\.\\w*$", "", names)

  treesL=vector(mode = "list", length = length(names))
  for(i in 1:length(alnfiles)){
    message(paste("Processing", alnfiles[i]))
    tree.res=estimatePhangornTree(alnfiles[i], treefile, submodel=submodel, type = type,
                                  format = format, k=k, ...)
    treesL[[i]]=tree.res$tree
  }
  fc=file(output.file, "wt")
  for(i in 1:length(alnfiles)){
    writeLines(paste(names[i], write.tree(treesL[[i]]), sep = "\t"), fc)
  }
  close(fc)
}
