#' Plot the provided tree, (optionally) rerooted, with specified branches highlighted
#'  
#' @param tree. A tree object.
#' @param outgroup. A vector of species to use to root the tree. If not provided, the tree remains unrooted.
#' @param hlspecies. A vector of species whose terminal branches to highlight, or a vector of branch numbers within the tree.
#' @param hlcols. Colors to use in highlighting the branches. If not specified, will use default R colors.
#' @param main. Main text for plot.
#' @return A plot of the the (optionally rerooted) tree, with branches highlighted.

plotTreeHighlightBranches <- function(tree, outgroup=NULL, hlspecies, hlcols=NULL, main=""){
  if (is.null(hlcols)) {
    hlcols <- c(2:(length(hlspecies)+1))
  }
  if (length(hlcols) < length(hlspecies)) {
    hlcols <- rep_len(hlcols, length(hlspecies))
  }
  if (!is.null(outgroup)) {
    outgroup <- outgroup[outgroup %in% tree$tip.label]
    if (length(outgroup) > 0) {
      #root the tree
      if (is.numeric(hlspecies)) {
        #track the branches
        dummytree <- tree
        dummytree$edge.length <- c(rep(1,nrow(dummytree$edge)))
        for (i in c(1:length(hlspecies))) {
          dummytree$edge.length[hlspecies] <- i+1
        }
        dummyrooted <- root(dummytree, outgroup)
      }
      rooted <- root(tree, outgroup)
    } else {
      print("No members of requested outgroup found in tree; keeping unrooted.")
      rooted <- tree
      outgroup <- NULL
    }
  } else {
    rooted <- tree
  }
  colMaster <- c(rep("black",nrow(rooted$edge)))
  if (is.numeric(hlspecies)) {
    if (!is.null(outgroup)) {
      hlcols <- c("black",hlcols)
      colMaster <- hlcols[dummyrooted$edge.length]
    } else {
      for (i in 1:length(hlspecies)) {
        colMaster[hlspecies[i]] <- hlcols[i]
      }
    }
  } else {
    wspmr <- rooted$tip.label[rooted$edge[,2]] #order of tips in edges
    for (i in 1:length(hlspecies)) {
      colMaster[which(wspmr==hlspecies[i])] <- hlcols[i]
    }
  }
  termedge <- order(rooted$edge[,2])[1:length(rooted$tip.label)] #which edge corresponds to each terminal branch
  colMasterTip <- colMaster[termedge]
  #Make branches of length 0 just *slightly* larger values to visualize tree
  rooted2=rooted
  mm=min(rooted2$edge.length[rooted2$edge.length>0])
  rooted2$edge.length[rooted2$edge.length==0]=max(0.02,mm/20)
  plotobj = plot.phylo(rooted2, main = main, edge.color=colMaster,
                       tip.color = colMasterTip, edge.width = 2, cex=0.8)
  return(plotobj)
}