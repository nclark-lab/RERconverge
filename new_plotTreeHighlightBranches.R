require(ggtree)
require(ape)

#' @param hlspecies the species to highlight
#' @param hlcols the colors to highlight those species with
#' @param main the name to give the plot
plotTreeHighlightBranches = function(tree, outgroup=NULL, hlspecies, hlcols=NULL, main="")
{
  if (!is.null(outgroup)) {
    outgroup <- outgroup[outgroup %in% tree$tip.label]
    if (length(outgroup) > 0) {
      #root the tree
      rooted <- root(tree, outgroup)
    } else {
      print("No members of requested outgroup found in tree; keeping unrooted.")
      rooted <- tree
      outgroup <- NULL
    }
  } else {
    rooted <- tree
  }
  
  #if hlcols is null make branches blue
  if (is.null(hlcols)) {
    hlcols <- rep_len("#0000ff", length(hlspecies)) #if there are fewer colors than species to highlight, repeat colors
  }else if (length(hlcols) < length(hlspecies)) {
    hlcols <- rep_len(hlcols, length(hlspecies)) 
  }
  
  #create hlspecies_named
  hlspecies_named <- vector(mode="character")
  if(is.numeric(hlspecies)){
    for(i in 1: length(hlspecies))
    {
      hlspecies_named[i] <- tree$tip.label[hlspecies[i]]
    }
  }else
    hlspecies_named <- hlspecies
  
  #Make branches of length 0 just *slightly* larger values to visualize tree
  rooted2 <- rooted
  mm <- min(rooted2$edge.length[rooted2$edge.length>0])
  rooted2$edge.length[rooted2$edge.length==0] <- max(0.02,mm/20)
  
 
  #vector of tip label colors
  tipCols <- vector(mode = "character")
  x <- 1
  for(i in 1: length(tree$tip.label))
  {
    if(tree$tip.label[i] %in% hlspecies_named)
    {
      tipCols[i] <- hlcols[x]
      x <- x+1
    }
    else tipCols[i] <- "black"
  }
  
  #number of labels
  nlabel <- rooted2$Nnode + length(rooted2$tip.label)
  
  edgeCols <- vector(mode="character", length=nlabel)
  x <- 1
  for(i in 1: nlabel)
  {
    if(i < length(rooted2$tip.label))
    {
      if(rooted2$tip.label[i] %in% hlspecies_named)
      {
        edgeCols[i] <- hlcols[x]
        x <- x+1
      }else
        edgeCols[i] <- "Black"
    }else
      edgeCols[i] <- "Black"
  }
  
  plotobj = ggtree(rooted2, color = edgeCols)
  plotobj = plotobj + geom_tiplab(color= tipCols, geom="text", cex = 3) + labs(title = main)
  
  return(plotobj)
}