#' RERconverge
#'
#'
#' @docType package
#' @author
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib RERconverge
#' @name RERconverge


require(ape)
require(phytools)
require(compiler)
require(plotrix)
require(Rcpp)
require(RcppArmadillo)
require(phangorn)
require(weights)
require("castor")
require(FSA)
require(Matrix)
require(data.table)

#' reads trees from a 2 column , tab seperated, file
#' The first columns is the gene name and the second column is the corresponding tree in parenthetic format known as the Newick or New Hampshire format

#' @param file The path to the tree file
#' @param  max.read This function takes a while for a whole genome, so max.read is useful for testing
#' @param  masterTree (optional) User can specify a master tree; only the topology will be used, and branch lengths will be inferred from gene trees.
#' @param  masterTree (optional) User can specify a master tree. Recommended only when
#' the number of available gene trees with all species is small.
#' @param  minTreesAll The minimum number of trees with all species present in order to estimate
#' master tree edge lengths (default 20).
#' @param reestimateBranches Boolean indicating whether to re-estimate branch lengths if master tree topology is included (default FALSE)
#' @param minSpecs the minimum number of species that needs to be present in a gene tree to be included in calculating master tree
#' @return A trees object of class "treeObj"
#' @export
readTrees=function(file, max.read=NA, masterTree=NULL, minTreesAll=20, reestimateBranches=F, minSpecs=NULL){
  tmp=scan(file, sep="\t", what="character", quiet = T)
  message(paste0("Read ",length(tmp)/2, " items", collapse=""))
  trees=vector(mode = "list", length = min(length(tmp)/2,max.read, na.rm = T))
  treenames=character()
  maxsp=0; # maximum number of species
  allnames=NA # unique tip labels in gene trees



  #create trees object, get species names and max number of species
  for ( i in 1:min(length(tmp),max.read*2, na.rm = T)){
    if (i %% 2==1){
      treenames=c(treenames, tmp[i])
    }
    else{
      trees[[i/2]]=tryCatch(unroot(read.tree(text=tmp[i])),
                            error = function(e) {
                              message('Cannot parse tree for the following gene: ',treenames[i/2]);
                              stop()
                            })

      #reduce to species present in master tree
      if (!is.null(masterTree)) {
        trees[[i/2]] = pruneTree(trees[[i/2]],intersect(trees[[i/2]]$tip.label,masterTree$tip.label))
      }

      # #check if it has more species
      # if(length(trees[[i/2]]$tip.label)>maxsp){
      #   maxsp=length(trees[[i/2]]$tip.label)
      #   allnames=trees[[i/2]]$tip.label
      # }

      #check if it has new species
      if (sum(trees[[i/2]]$tip.label %in% allnames == F) > 0) {
        allnames = unique(c(allnames,trees[[i/2]]$tip.label))
        maxsp = length(allnames) - 1

      }
      #if(length(trees[[i/2]]$tip.label)>maxsp){
      #  maxsp=length(trees[[i/2]]$tip.label)
      #  allnames=trees[[i/2]]$tip.label
      #}

    }

  }


  allnames = allnames[!is.na(allnames)]
  names(trees)=treenames
  treesObj=vector(mode = "list")
  treesObj$trees=trees
  treesObj$numTrees=length(trees)
  treesObj$maxSp=maxsp

  message(paste("max is", maxsp))



  ### report is a binary matrix showing the species membership of each tree
  report=matrix(nrow=treesObj$numTrees, ncol=maxsp)
  colnames(report)=allnames

  rownames(report)=treenames
  for ( i in 1:nrow(report)){
    ii=match(allnames, trees[[i]]$tip.label)
    report[i,]=1-is.na(ii)

  }
  treesObj$report=report

  ############ This line finds indices of trees that have the complete species
  ii=which(rowSums(report)==maxsp)



  ######################################################################
  if(length(ii)==0 & is.null(masterTree)){
    stop("no tree has all species - you must supply a master tree")
  }
  ######################################################################


  #Create a master tree with no edge lengths
  if (is.null(masterTree)) {
    master=trees[[ii[1]]]
    master$edge.length[]=1
    treesObj$masterTree=master
  } else {

    master=pruneTree(masterTree, intersect(masterTree$tip.label,allnames))
    #prune tree to just the species names in the largest gene tree
    master$edge.length[]=1

    master=unroot(pruneTree(masterTree, intersect(masterTree$tip.label,allnames)))
    #prune tree to just the species names in the gene trees
    #master$edge.length[]=1

    treesObj$masterTree=master
  }


  treesObj$masterTree=rotateConstr(treesObj$masterTree, sort(treesObj$masterTree$tip.label))
  #this gets the abolute alphabetically constrained order when all branches
  #are present
  tiporder=treeTraverse(treesObj$masterTree)

  #treesObj$masterTree=CanonicalForm(treesObj$masterTree)
  message("Rotating trees")

  for ( i in 1:treesObj$numTrees){

    treesObj$trees[[i]]=rotateConstr(treesObj$trees[[i]], tiporder)

  }


  ap=allPathsTrackBranches(master)
  treesObj$ap=ap
  matAnc=(ap$matIndex>0)+1-1
  matAnc[is.na(matAnc)]=0



  paths=matrix(nrow=treesObj$numTrees, ncol=length(ap$dist))
  for( i in 1:treesObj$numTrees){
    #Make paths all NA if tree topology is discordant
    paths[i,]=tryCatch(allPathMasterRelativeTrackBranches(treesObj$trees[[i]], master, ap), error=function(err) NA)
    #calls matchAllNodes -> matchNodesInject
  }
  paths=paths+min(paths[paths>0], na.rm=T)
  treesObj$paths=paths
  treesObj$matAnc=matAnc
  treesObj$matIndex=ap$matIndex
  treesObj$lengths=unlist(lapply(treesObj$trees, function(x){sqrt(sum(x$edge.length^2))}))

  #require all species and tree compatibility
  #ii=which(rowSums(report)==maxsp)
  ii=intersect(which(rowSums(report)==maxsp),which(is.na(paths[,1])==FALSE))



  #if masterTree is provided by user, must use minSpecs<maxsp
  #if no user supplied tree and not minSpec, calculate branch lengths from trees with all species
  #if minSpecs<maxsp, calculate branch lengths from trees with minSpecs species
  if(is.null(minSpecs)){
    #if minimum species not specified,
    #minimum is all species
    minSpecs=maxsp
  }

  if(!is.null(masterTree) && !reestimateBranches){
    message("Using user-specified master tree")
  }


  if(minSpecs==maxsp){ #if we're using all species
    if (is.null(masterTree)) { #and if the user did not specify a master tree
      if(length(ii)>=minTreesAll){
        message (paste0("estimating master tree branch lengths from ", length(ii), " genes"))
        tmp=lapply( treesObj$trees[ii], function(x){x$edge.length})

        allEdge=matrix(unlist(tmp), ncol=2*maxsp-3, byrow = T)
        allEdge=scaleMat(allEdge)
        allEdgeM=apply(allEdge,2,mean)
        treesObj$masterTree$edge.length=allEdgeM
      }else {
        message("Not enough genes with all species present: master tree has no edge.lengths")
      }
    }else{
      message("Must specify minSpecs when supplying a master tree: master tree has no edge.lengths")
    }
  }else{ #if we are not using all species
    #estimating from trees with minimum number of species
    treeinds=which(rowSums(report)>=minSpecs) #which trees have the minimum species
    message (paste0("estimating master tree branch lengths from ", length(treeinds), " genes"))


    if(length(treeinds)>=minTreesAll){
      pathstouse=treesObj$paths[treeinds,] #get paths for those trees

      colnames(pathstouse) = ap$destinNode
      colBranch = vector("integer",0)
      unq.colnames = unique(colnames(pathstouse))

      for (i in 1:length(unq.colnames)){
        ind.cols = which(colnames(pathstouse) == unq.colnames[i])
        colBranch = c(colBranch,ind.cols[1])
      }

      allEdge = pathstouse[,colBranch]
      allEdgeScaled = allEdge
      for (i in 1:nrow(allEdgeScaled)){
        allEdgeScaled[i,] = scaleDistNa(allEdgeScaled[i,])
      }
      colnames(allEdgeScaled) = unq.colnames

      edgelengths = vector("double", ncol(allEdgeScaled))

      edge.master = treesObj$masterTree$edge

      for (i in 1:nrow(edge.master)){
        destinNode.i = edge.master[i,2]
        col.Node.i = allEdgeScaled[,as.character(destinNode.i)]
        edgelengths[i] = mean(na.omit(col.Node.i))
      }

      treesObj$masterTree$edge.length = edgelengths
    }else{
      message("Not enough genes with minSpecs species present: master tree has no edge.lengths")
    }
  }

  message("Naming columns of paths matrix")
  colnames(treesObj$paths)=namePathsWSpecies(treesObj$masterTree)
  class(treesObj)=append(class(treesObj), "treesObj")
  treesObj
}

#' @keywords  internal
allPathMasterRelativeTrackBranches=function(tree, masterTree, masterTreePaths=NULL){
  if(! is.list(masterTreePaths)){
    masterTreePaths=allPathsTrackBranches(masterTree)
  }

  treePaths=allPaths(tree)
  map=matchAllNodes(tree,masterTree)

  #remap the nodes
  treePaths$nodeId[,1]=map[treePaths$nodeId[,1],2 ]
  treePaths$nodeId[,2]=map[treePaths$nodeId[,2],2 ]


  ii=masterTreePaths$matIndex[(treePaths$nodeId[,2]-1)*nrow(masterTreePaths$matIndex)+treePaths$nodeId[,1]]

  vals=double(length(masterTreePaths$dist))
  vals[]=NA
  vals[ii]=treePaths$dist
  vals
}

#' @keywords  internal
allPathsTrackBranches=function(tree){
  dd=dist.nodes(tree) #### pairwise distances between nodes in the tree
  allD=double()
  nn=matrix(nrow=0, ncol=2)
  nA=length(tree$tip.label)+tree$Nnode ######### Total number of nodes in the tree (internal + tips)
  matIndex=matrix(nrow=nA, ncol=nA)
  index=1

  destinNode = vector("integer", 0)
  ancNode = vector("integer",0)

  for ( i in 1:nA){
    ia=getAncestors(tree,i) #### Getting the ancestors of each node in the tree

    destinNode = c(destinNode, rep(i, length(ia)))
    ancNode = c(ancNode, ia)

    if(length(ia)>0){
      allD=c(allD, dd[i, ia])
      nn=rbind(nn,cbind(rep(i, length(ia)), ia))
      for (j in ia){
        matIndex[i,j]=index
        index=index+1
      }
    }
  }
  return(list(dist=allD, nodeId=nn, matIndex=matIndex, destinNode=destinNode, ancNode=ancNode))
}

#' @keywords  internal
scaleDistNa=function(x){
  x/sqrt(sum(x[!is.na(x)]^2))
}







computeWeightsAllVar=function (mat, nv=NULL, transform="none",plot = T, predicted=T){

  if(is.null(nv)){
    nv=apply(mat, 2, mean,na.rm=T, trim=0.05)

  }
  transform=match.arg(transform, choices = c("none", "sqrt", "log"))
  if(transform=="log"){
    offset=0;
    if(min(mat, na.rm=T)<1e-8){
      offset=min(mat[mat>1e-8])
    }
    mat=log(mat+offset)
    nv=log(nv+offset)
  }
  if (transform=="sqrt"){
    mat=sqrt(mat)
    nv=sqrt(nv)
  }

  matsub=mat


  matr=naresidCPP(matsub, model.matrix(~1+nv))
  matpred=fastLmPredictedMat(matsub, model.matrix(~1+nv))



  mml=as.vector(matsub)
  varl=as.vector(log(matr^2))
  ii=which(!is.na(mml))
  mml=mml[ii]
  varl=varl[ii]
  set.seed(123)
  iis=sample(length(mml), min(500000, length(mml)))
  mml=mml[iis]
  varl=varl[iis]

  l = lowess(mml,varl, f=0.7, iter = 2)

  f = approxfun(l, rule = 2)
  if (plot) {
    par(mfrow=c(1,2), omi=c(1,0,0,0))
    nbreaks=20
    qq=quantile(mml,seq(0,nbreaks,1)/nbreaks)
    qqdiff=diff(qq)
    breaks=qq[1:nbreaks]+qqdiff/2
    rr=quantile(mml, c(0.0001, 0.99))
    #breaks=round(breaks,3)
    breaks=unique(round(breaks,3)) #forces unique breaks
    nbreaks=length(breaks)
    cutres<-cut(mml,breaks = breaks)

    cutres_tt=table(cutres)

    boxplot((varl)~ cutres, xlab = "", ylab = "log var", outline=F,  log="", las=2)
    title("Before")


    xx=(qq[1:nbreaks]+breaks)/2

    lines(1:length(xx), (f(qq[1:nbreaks])), lwd = 2, col = 2)
  }
  wr=1/exp(f(mml))

  if(!predicted){
    weights=(matrix(1/exp(f(mat)), nrow = nrow(mat)))
  }
  else{
    weights=(matrix(1/exp(f(matpred)), nrow = nrow(mat)))
  }
  if(plot){
    matr=naresidCPP(matsub, model.matrix(~1+nv), weights)
    varl=(as.vector(log(matr^2))[ii])[iis]
    boxplot((varl)~ cutres, ylab = "log var", outline=F,  log="", main="After", las=2)
    abline(h=0, col="blue3", lwd=2)
    mtext(side = 1, text="bins", outer = T, line = 2)
  }
  weights
}

#' @keywords  internal
naresidCPP=function(data, mod, weights=NULL){
  if(is.null(weights)){
    out=fastLmResidMat(data, mod)
  }
  else{
    out=fastLmResidMatWeighted(data,mod, weights)
  }
  rownames(out)=rownames(data)
  colnames(out)=colnames(data)
  out
}

#' Provides names for paths/RERs representing terminal branches for plotting
#' Originally an internal function but necessary for the vignette/walk-through
#' @param  masterTree The master tree used for analysis
#' @return  Names corresponding to the paths/RERs for terminal branches
#' @export
namePathsWSpecies=function(masterTree){
  mat=transformMat(masterTree)
  n=length(masterTree$tip.label)
  #these are the tip edges in the master tree
  iim=match(1:n, masterTree$edge[,2])
  #each column in the mat is composed of at most one tip edge
  tip.edge=apply(mat[iim,],2,function(x){if(max(x)>0){which(x==1)} else{NA}})
  return(masterTree$tip[tip.edge])

}
#' @keywords  internal
scaleMat=function(mat){t(apply(mat,1,scaleDist))}



#' @keywords  internal
scaleDist=function(x){
  x/sqrt(sum(x^2))
}

#' @keywords  internal
allPathMasterRelative=function(tree, masterTree, masterTreePaths=NULL){
  if(! is.list(masterTreePaths)){
    masterTreePaths=allPaths(masterTree)
  }

  treePaths=allPaths(tree)
  map=matchAllNodes(tree,masterTree)

  #remap the nodes
  treePaths$nodeId[,1]=map[treePaths$nodeId[,1],2 ]
  treePaths$nodeId[,2]=map[treePaths$nodeId[,2],2 ]


  ii=masterTreePaths$matIndex[(treePaths$nodeId[,2]-1)*nrow(masterTreePaths$matIndex)+treePaths$nodeId[,1]]

  vals=double(length(masterTreePaths$dist))
  vals[]=NA
  vals[ii]=treePaths$dist
  vals
}





#' @keywords  internal
matchAllNodes=function(tree1, tree2){
  map=matchNodesInject(tree1,tree2)
  map=map[order(map[,1]),]
  map
}


#' @keywords  internal
matchNodesInject=function (tr1, tr2){
  if(length(tmpsp<-setdiff(tr1$tip.label, tr2$tip.label))>0){
    #stop(paste(paste(tmpsp, ","), "in tree1 do not exist in tree2"))
    stop(c("The following species in tree1 do not exist in tree2: ",paste(tmpsp, ", ")))
  }
  commontiplabels <- intersect(tr1$tip,tr2$tip)
  if(RF.dist(pruneTree(tr1,commontiplabels),pruneTree(tr2,commontiplabels))>0){
    stop("Discordant tree topology detected - gene/trait tree and treesObj$masterTree have irreconcilable topologies")
  }
  #if(RF.dist(tr1,tr2)>0){
  #  stop("Discordant tree topology detected - trait tree and treesObj$masterTree have irreconcilable topologies")
  #}

  toRm=setdiff(tr2$tip.label, tr1$tip.label)
  desc.tr1 <- lapply(1:tr1$Nnode + length(tr1$tip), function(x) extract.clade(tr1,
                                                                              x)$tip.label)
  names(desc.tr1) <- 1:tr1$Nnode + length(tr1$tip)
  desc.tr2 <- lapply(1:tr2$Nnode + length(tr2$tip), function(x) extract.clade(tr2,
                                                                              x)$tip.label)
  names(desc.tr2) <- 1:tr2$Nnode + length(tr2$tip)
  Nodes <- matrix(NA, length(desc.tr1), 2, dimnames = list(NULL,
                                                           c("tr1", "tr2")))
  for (i in 1:length(desc.tr1)) {
    Nodes[i, 1] <- as.numeric(names(desc.tr1)[i])
    for (j in 1:length(desc.tr2)) if (all(desc.tr1[[i]] %in%
                                          desc.tr2[[j]]))
      Nodes[i, 2] <- as.numeric(names(desc.tr2)[j])
  }

  iim=match(tr1$tip.label, tr2$tip.label)
  Nodes=rbind(cbind(1:length(tr1$tip.label),iim),Nodes)
  if(any(table(Nodes[,2])>1)){
    stop("Incorrect pseudorooting detected - use fixPseudoroot() function to correct trait tree topology")
  }

  Nodes
}



#' @keywords  internal
allPaths=function(tree, categorical = F){
  if (!categorical){
    dd=dist.nodes(tree)
  }
  allD=double()
  nn=matrix(nrow=0, ncol=2)
  nA=length(tree$tip.label)+tree$Nnode
  matIndex=matrix(nrow=nA, ncol=nA)
  index=1
  for ( i in 1:nA){
    ia=getAncestors(tree,i)
    if(length(ia)>0){
      if(categorical) {
        # add the state of node i to allD length(ia) times
        x = which(tree$edge[,2] == i)
        state = tree$edge.length[x]
        allD = c(allD, rep(state, length(ia)))
      }
      else {
        allD=c(allD, dd[i, ia])
      }
      nn=rbind(nn,cbind(rep(i, length(ia)), ia))
      for (j in ia){
        matIndex[i,j]=index
        index=index+1
      }
    }
  }
  return(list(dist=allD, nodeId=nn, matIndex=matIndex))
}

#' @keywords  internal
getAncestors=function(tree, nodeN){
  if(is.character(nodeN)){
    nodeN=which(tree$tip.label==nodeN)
  }
  im=which(tree$edge[,2]==nodeN)
  if(length(im)==0){
    return()
  }
  else{
    anc=tree$edge[im,1]
    return(c(anc, getAncestors(tree, anc)))
  }

}




#' @keywords  internal
treeTraverse=function(tree, node=NULL){
  if(is.null(node)){
    rt=getRoot(tree)
    ic=getChildren(tree,rt)
    return(c(treeTraverse(tree, ic[1]), treeTraverse(tree, ic[2])))

  }
  else{
    if (node<=length(tree$tip)){
      return(tree$tip[node])
    }
    else{
      ic=getChildren(tree,node)
      return(c(treeTraverse(tree, ic[1]), treeTraverse(tree, ic[2])))

    }
  }
}
#' @keywords  internal
getRoot = function(phy) phy$edge[, 1][!match(phy$edge[, 1], phy$edge[, 2], 0)][1]
#' @keywords  internal
getChildren=function(tree, nodeN){
  tree$edge[tree$edge[,1]==nodeN,2]
}



#'Computes the association statistics between RER from \code{\link{getAllResiduals}} and a phenotype paths vector for a binary phenotype made with \code{\link{tree2Paths}}
#' @param RERmat RER matrix returned by \code{\link{getAllResiduals}}
#' @param charP phenotype vector returned by \code{\link{tree2Paths}} or \code{\link{char2Paths}}
#' @param method Method used to compute correlations. Accepts the same arguments as \code{\link{cor}}. Set to "auto" to select automatically based on the number of unique values in charP. This will also auto set the winsorization for Pearson correlation. Set winsorizetrait=some number and winsorizeRER=some number to override
#' @param min.sp Minimum number of species that must be present for a gene
#' @param min.pos Minimum number of species that must be present in the foreground (non-zero phenotype values)
#' @param weighted perform weighted correlation. This option turns on weighted correlation that uses the weights computed by \code{\link{foreground2Tree}(wholeClade=T)}. This setting will treat each clade a single observation for the purpose of p-value estimation. The function will guess automatically if the charP vector is of "weighted" type and there should be not need to set this parameter.
#' @export
correlateWithBinaryPhenotype=function(RERmat,charP, min.sp=10, min.pos=2, weighted="auto"){
  if(weighted=="auto"){
    if (any(charP>0&charP<1, na.rm=TRUE)){
      message("Fractional values detected, will use weighted correlation mode")
      weighted=T
    }
    else{
      weighted=F
    }
  }
  getAllCor(RERmat, charP, min.sp, min.pos, method = "k", weighted=weighted)

}


#'Computes the association statistics between RER from \code{\link{getAllResiduals}} and a phenotype paths vector for a continuous phenotype made with \code{\link{char2Paths}}
#' @param RERmat RER matrix returned by \code{\link{getAllResiduals}}
#' @param charP phenotype vector returned by \code{\link{tree2Paths}} or \code{\link{char2Paths}}
#' @param method Method used to compute correlations. Accepts the same arguments as \code{\link{cor}}. Set to "auto" to select automatically based on the number of unique values in charP. This will also auto set the winsorization for Pearson correlation. Set winsorizetrait=some number and winsorizeRER=some number to override
#' @param min.sp Minimum number of species that must be present for a gene
#' @param winsorizeRER Winsorize RER values before computing Pearson correlation. winsorizeRER=3 will set the 3 most extreme values at each end of each RER row to the the value closest to 0.
#' @param winsorizetrait Winsorize trait values before computing Pearson correlation. winsorizetrait=3 will set the 3 most extreme values of the trait values to the value closest to 0.
#' @export

correlateWithContinuousPhenotype=function(RERmat,charP, min.sp=10,  winsorizeRER=3, winsorizetrait=3){
  getAllCor(RERmat, charP, min.sp, min.pos=0, method = "p", winsorizeRER = winsorizeRER, winsorizetrait = winsorizetrait)
}


#' Computes the association statistics between RER from \code{\link{getAllResiduals}} and a phenotype paths vector for a categorical phenotype made with \code{\link{char2PathsCategorical}}
#'@param RERmat RER matrix returned by \code{\link{getAllResiduals}}
#'@param charP phenotype vector returned by \code{\link{tree2Paths}} or \code{\link{char2Paths}}
#'@param method Method used to compute correlations. Use "kw" to use Kruskil Wallis. Use "aov" to use ANOVA. It must be one of the strings "kw" or "aov".
#'@param min.sp Minimum number of species that must be present for a gene
#'@param min.pos Minimum number of species that must be present in a category
#'@return A list containing a list object with correlation values, p-values, and the number of data points used for each tree and a list of list objects for each pairwise test with correlation values and p-values.
#'@export
correlateWithCategoricalPhenotype = function(RERmat,charP, min.sp = 10, min.pos = 2, method = "kw"){
  if(!(method %in% c("kw", "aov"))) {
    warning("Invalid method. The method must be kw or aov")
  }
  getAllCor(RERmat, charP, min.sp, min.pos, method = method)
}

#' A sped up version of the Kruskal Wallis/Dunn Test 
#' @keywords internal
kwdunn.test <- function(x,g, ncategories){
  ntests <- ncategories*(ncategories -1)/2 # number of pairwise tests
  # set up a data matrix
  N <- length(x)
  glevels <- as.integer(levels(g)) # g must be an integer converted to a factor!
  Data <- matrix(NA, length(x), 3)
  Data[, 1] <- x
  Data[, 2] <- g  
  # when g gets converted to numeric in the Data matrix, it is converted to CONSECUTIVE integers in the order of the factors
  # e.g. 1, 2, 4 --> 1, 2, 3
  
  # use frank (fast rank) instead of base rank function
  # REQUIRES THE PACKAGE data.table TO BE ATTACHED!!!
  Data[, 3] <- frank(Data[, 1], ties.method = "average", na.last = NA)
  
  # calculate the ties adjustment term that is shared between kwallis and dunn test
  # define a function to find the tied ranks
  tiedranks <- function(ranks) {
    ranks <- sort(ranks)
    ties <- c()
    for (i in 2:length(ranks)) {
      if (ranks[i - 1] == ranks[i]) {
        if (length(ties) > 0) {
          if (ranks[i - 1] != tail(ties, n = 1)) {
            ties <- c(ties, ranks[i - 1])
          }
        }
        else {
          ties <- c(ranks[i - 1])
        }
      }
    }
    return(ties)
  }
  
  # calculate the ties adjustment sum term (same between KW and Dunn)
  k <- length(unique(Data[, 2]))
  ranks <- Data[, 3]
  ties <- tiedranks(ranks)
  r <- length(ties)
  tiesadjsum <- 0
  if (r > 0) {
    for (s in 1:r) {
      tau <- sum(ranks == ties[s])
      tiesadjsum <- tiesadjsum + (tau^{3} - tau)
    }
  }
  
  # pre-calculate indices/sums/stuff to reduce the number of times it's calculated
  groupinds <- lapply(1:k, function(i){Data[, 2] == i}) # rows in Data corresponding to group i (as TRUE/FALSE vector)
  groupranks <- lapply(groupinds, function(i){Data[, 3][i]}) # ranks in each group
  groupranksums <- unlist(lapply(groupranks, function(i){sum(i)})) # sum of ranks in each group
  groupsizes <- unlist(lapply(groupinds, function(i){sum(i)})) # number of observations in each group
  
  # calculate the H statistic and p-value for the KW test
  tiesadj <- 1 - (tiesadjsum/((N^3) - N))
  ranksum <- sum((groupranksums^2)/groupsizes) # use matrix operations in place of for loops
  H <- ((12/(N * (N + 1))) * ranksum - 3 * (N + 1))/tiesadj
  df <- k - 1
  p <- pchisq(H, k - 1, lower.tail = FALSE)
  
  # Dunn test: calculate the Z statistic for each pairwise test 
  m <- k * (k - 1)/2
  Z <- rep(NA, ntests)
  tiesadj <- tiesadjsum/(12 * (N - 1))
  
  # loop through each pairwise comparison
  index <- 1
  for (i in 2:k) {
    for (j in 1:(i - 1)) {
      # make a pairwise test name
      index <- ((glevels[i]-1) * (glevels[i] - 2)/2) + glevels[j]
      
      # do calculation
      meanranki <- groupranksums[i]/groupsizes[i]
      meanrankj <- groupranksums[j]/groupsizes[j]
      z <- (meanrankj - meanranki)/sqrt(((N * (N + 1)/12) - tiesadj) * ((1/groupsizes[j]) + (1/groupsizes[i])))
      
      # add result to Z vector with the name
      Z[index] <- z
    }
  }
  P <- 2*pnorm(abs(Z), lower.tail = FALSE)
  
  # do bonferroni p value adjustment (for pairwise error rate?)
  P.adjust <- pmin(1, P * m)
  return(list(kw = list(H = H, p = p), dunn = list(Z = Z, P = P, P.adjust = P.adjust)))
}

#'Computes the association statistics between RER from \code{\link{getAllResiduals}} and a phenotype paths vector made with \code{\link{tree2Paths}} or \code{\link{char2Paths}}
#' @param RERmat RER matrix returned by \code{\link{getAllResiduals}}
#' @param charP phenotype vector returned by \code{\link{tree2Paths}} or \code{\link{char2Paths}}
#' @param method Method used to compute correlations. Accepts the same arguments as \code{\link{cor}}. Set to "auto" to select automatically based on the number of unique values in charP. This will also auto set the winsorization for Pearson correlation. Set winsorizetrait=some number and winsorizeRER=some number to override
#' @param min.sp Minimum number of species that must be present for a gene
#' @param min.pos Minimum number of species that must be present in the foreground (non-zero phenotype values)
#' @param winsorizeRER Winsorize RER values before computing Pearson correlation. winsorizeRER=3 will set the 3 most extreme RER values at each end of each row to the value closest to 0.
#' @param winsorizetrait Winsorize trait values before computing Pearson correlation. winsorizetrait=3 will set the 3 most extreme trait values at each end to the value closest to 0.
#' @param weighted perform weighted correlation. This option needs to be set if the clade weights computed in \code{\link{foreground2Tree}(wholeClade=T)} are to be used. This setting will treat the clade a single observation for the purpose of p-value estimation.
#' @note  winsorize is in terms of number of observations at each end, NOT quantiles
#' @return A list object with correlation values, p-values, and the number of data points used for each tree
#' @export
getAllCor=function(RERmat, charP, method="auto",min.sp=10, min.pos=2, winsorizeRER=NULL, winsorizetrait=NULL, weighted=F){
  RERna=(apply(is.na(RERmat),2,all))
  iicharPna=which(is.na(charP))
  if(!all(RERna[iicharPna])){
    warning("Species in phenotype vector are a subset of the those used for RER computation. For best results run getAllResiduals with the useSpecies")
  }
  if (method=="auto"){
    lu=length(unique(charP))
    if(lu==2){
      method="k"
      message("Setting method to Kendall")
    }
    else if (lu<=5){
      method="s"
      message("Setting method to Spearman")
    }
    else{
      method="p"
      message("Setting method to Pearson")
      if(is.null(winsorizeRER)){
        message("Setting winsorizeRER=3")
        winsorizeRER=3
      }
      if(is.null(winsorizetrait)){
        message("Setting winsorizetrait=3")
        winsorizetrait=3
      }
    }
  }
  win=function(x,w){
    xs=sort(x[!is.na(x)], decreasing = T)
    xmax=xs[w]
    xmin=xs[length(xs)-w+1]

    x[x>xmax]=xmax
    x[x<xmin]=xmin
    x
  }
  corout=matrix(nrow=nrow(RERmat), ncol=3)
  rownames(corout)=rownames(RERmat)

  ##############################################################################
  # make tables for each pairwise comparison & a list for the tables
  if (method == "aov" || method == "kw") {
    lu = length(unique(charP[!is.na(charP)]))
    n = choose(lu, 2)
    tables = lapply(1:n, matrix, data = NA, nrow = nrow(RERmat), 
                    ncol = 2, dimnames = list(rownames(RERmat), c("Rho", 
                                                                  "P")))
    if(method == "aov") {
      names(tables) = rep(NA, n)
    }
    else { # name the tables in the same order as Z is calculated
      for (i in 2:lu) {
        for (j in 1:(i - 1)) {
          index <- (i-1)*(i-2)/2 + j
          names(tables)[index] <- paste0(j, " - ", i)
        }
      }
    }
    
  }
  ##############################################################################

  colnames(corout)=c("Rho", "N", "P")

  for( i in 1:nrow(corout)){

    if(((nb<-sum(ii<-(!is.na(charP)&!is.na(RERmat[i,]))))>=min.sp)){
      ##########################################################################
      #check that there are >= min.pos species in each category
      if(method =="kw" || method =="aov") {
        counts = table(charP[ii])
        num_groups = length(counts) # k in eta^2 calculation for KW test
        if(num_groups < 2 || min(counts) < min.pos) {
          next
        }
      }
      ##########################################################################
      else if (method!="p"&&sum(charP[ii]!=0)<min.pos){
        next
      }
      if(!weighted){

        x=RERmat[i,]

        #winsorize
        indstouse=which(!is.na(x) & !is.na(charP))
        if(!is.null(winsorizeRER)){
          x=win(x[indstouse], winsorizeRER)
        }else{
          x=x[indstouse]
        }
        if(!is.null(winsorizetrait)){
          y=win(charP[indstouse], winsorizetrait)
        }else{
          y=charP[indstouse]
        }

        #################################################################
        if(method == "aov") {
          # ANOVA
          # make a data frame from x and y
          yfacts = as.factor(y) #aov requires categories to be factors
          df = data.frame(x,yfacts)
          colnames(df) = c("RER", "category")
          ares = aov(RER ~ category, data = df)
          # ares_Fval = summary(ares)[[1]][1,4]
          # calculate effect size:
          sumsq = summary(ares)[[1]][1,2]
          sumsqres = summary(ares)[[1]][2,2]
          # eta2
          effect_size = sumsq / (sumsq + sumsqres)
          ares_pval = summary(ares)[[1]][1,5]
          corout[i,1:3]=c(effect_size, nb, ares_pval)

          tukey = TukeyHSD(ares)

          # add names to tables that haven't been named yet
          groups = rownames(tukey[[1]])
          unnamedinds = which(is.na(names(tables)))
          if(length(unnamedinds > 0)) {
            # check for groups not in table already
            newnamesinds = which(is.na(match(groups, names(tables))))
            # if there are new groups add them to the next available positions
            if(length(newnamesinds) > 0) {
              names(tables)[unnamedinds][1:length(newnamesinds)] = groups[newnamesinds]
            }
          }

          # add data to the named tables
          for(k in 1:length(groups)) {
            name = groups[k]
            tables[[name]][i,"Rho"] = tukey[[1]][name,1]
            tables[[name]][i,"P"] = tukey[[1]][name,4]
          }

        } else if (method == "kw") {
          # Kruskal Wallis/Dunn test
          yfacts = factor(y)
          kres = kwdunn.test(x, yfacts, ncategories = lu)
          effect_size = kres$kw$H/(nb - 1)
          corout[i, 1:3] = c(effect_size, nb, kres$kw$p)
         
          for(k in 1:length(kres$dunn$Z)){ # length of kres$dunn$Z should be the same as length(tables) otherwise there's a problem
            tables[[k]][i, "Rho"] = kres$dunn$Z[k]
            tables[[k]][i, "P"] = kres$dunn$P.adjust[k]
          }
          
          # old code before speed up:
          # yfacts = as.factor(y)
          # df = data.frame(x,yfacts)
          # colnames(df) = c("RER", "category")
          # kres = kruskal.test(RER ~ category, data = df)
          # kres_Hval = kres$statistic
          # kres_pval = kres$p.value
          # # calculate effect size
          # # effect_size = (kres_Hval - num_groups + 1) / (nb - num_groups) # eta2: (H - k + 1) / (n - k)
          # effect_size = kres_Hval / (nb - 1) # epsilon squared
          # corout[i,1:3] = c(effect_size, nb, kres_pval)
          # 
          # # Dunn test
          # dunn = dunnTest(RER ~ category, data = df, method = "bonferroni") # do we want to use bonferroni?
          # 
          # # add new names to tables
          # groups = dunn$res$Comparison
          # unnamedinds = which(is.na(names(tables)))
          # if(length(unnamedinds > 0)) {
          #   # check for groups not in table already
          #   newnamesinds = which(is.na(match(groups, names(tables))))
          #   # if there are new groups add them to the next available positions
          #   if(length(newnamesinds) > 0) {
          #     names(tables)[unnamedinds][1:length(newnamesinds)] = groups[newnamesinds]
          #   }
          # }
          # # add data to the tables
          # for(k in 1:length(groups)) {
          #   name = groups[k]
          #   tables[[name]][i,"Rho"] = dunn$res$Z[k]
          #   tables[[name]][i,"P"] = dunn$res$P.adj[k]
          # }
        }
        else {
          cres=cor.test(x, y, method=method, exact=F)
          corout[i,1:3]=c(cres$estimate, nb, cres$p.value)
        }
      }
      else{
        charPb=(charP[ii]>0)+1-1

        weights=charP[ii]
        weights[weights==0]=1

        cres=wtd.cor(RERmat[i,ii], charPb, weight = weights, mean1 = F)
        corout[i, 1:3]=c(cres[1], nb, cres[4])
      }
    }
    else{
      #show(i)
      #show(c(nb, charP[ii]))
    }

  }

  corout=as.data.frame(corout)
  corout$p.adj=p.adjust(corout$P, method="BH")
  if (method == "aov" || method == "kw") {
    for(i in 1:length(tables)) {
      tables[[i]] = as.data.frame(tables[[i]])
      # add adjusted p-values
      tables[[i]]$p.adj = p.adjust(tables[[i]]$P, method = "BH")
    }
    # return corout and tables
    return(list(corout,tables))
  } else {corout}
}

#'Computes the association statistics between RER from \code{\link{getAllResiduals}} and a phenotype vector for phenotype values at the tips of the tree
#' @param RERmat RER matrix returned by \code{\link{getAllResiduals}}
#' @param method Method used to compute correlations. Accepts the same arguments as \code{\link{cor}}. Set to "auto" to select automatically based on the number of unique values in charP. This will also auto set the winsorization for Pearson correlation. Set winsorizetrait=some number and winsorizeRER=some number to override
#' @param min.sp Minimum number of species that must be present for a gene
#' @param min.pos Minimum number of species that must be present in the foreground (non-zero phenotype values)
#' @param winsorizeRER Winsorize RER values before computing Pearson correlation. winsorizeRER=3 will set the 3 most extreme RER values at each end of each row to the value closest to 0.
#' @param winsorizetrait Winsorize trait values before computing Pearson correlation. winsorizetrait=3 will set the 3 most extreme trait values at each end to the value closest to 0.
#' @note  winsorize is in terms of number of observations at each end, NOT quantiles
#' @return A list object with correlation values, p-values, and the number of data points used for each tree
#' @export
getAllCorExtantOnly <- function (RERmat, phenvals, method = "auto",
                                 min.sp = 10, min.pos = 2, winsorizeRER = NULL,
                                 winsorizetrait = NULL)
{
  # if method is auto, set the method
  if (method == "auto") {
    # lu = length(unique(charP))
    lu = length(unique(phenvals))
    if (lu == 2) {
      method = "k"
      message("Setting method to Kendall")
    }
    else if (lu <= 5) {
      method = "s"
      message("Setting method to Spearman")
    }
    else {
      method = "p"
      message("Setting method to Pearson")
      if (is.null(winsorizeRER)) {
        message("Setting winsorizeRER=3")
        winsorizeRER = 3
      }
      if (is.null(winsorizetrait)) {
        message("Setting winsorizetrait=3")
        winsorizetrait = 3
      }
    }
  }

  # define the function, win
  win = function(x, w) {
    xs = sort(x[!is.na(x)], decreasing = T)
    xmax = xs[w]
    xmin = xs[length(xs) - w + 1]
    x[x > xmax] = xmax
    x[x < xmin] = xmin
    x
  }

  # make the matrix to store the results
  corout = matrix(nrow = nrow(RERmat), ncol = 3)
  rownames(corout) = rownames(RERmat)

  # generate tables for the pairwise tests if trait is categorical
  if (method == "aov" || method == "kw") {
    # convert phenvals to integers to work with kwdunn.test function
    intlabels = map_to_state_space(phenvals)
    phenv = intlabels$mapped_states
    names(phenv) = names(phenvals)
    phenvals = phenv
    
    lu = length(unique(phenvals[!is.na(phenvals)]))
    n = choose(lu, 2)
    tables = lapply(1:n, matrix, data = NA, nrow = nrow(RERmat), 
                    ncol = 2, dimnames = list(rownames(RERmat), c("Rho", 
                                                                  "P")))
    if(method == "aov") {
      names(tables) = rep(NA, n)
    }
    else { # name the tables in the same order as Z is calculated
      for (i in 2:lu) {
        for (j in 1:(i - 1)) {
          index <- (i-1)*(i-2)/2 + j
          names(tables)[index] <- paste0(j, " - ", i)
        }
      }
    }
  }

  # name the columns of corout
  colnames(corout) = c("Rho", "N", "P")

  # for each gene in the analysis...
  for (i in 1:nrow(corout)) {
    rer = RERmat[i,]
    # get rid of NA values
    rer = rer[!is.na(rer)]
    # get rid of unnamed values (internal nodes)
    rer = rer[!is.na(names(rer))]

    # get the groups
    phens = phenvals

    # find which species are in phens and rers
    keep = intersect(names(phens), names(rer))
    # remove the rer values for species not in phens
    rer = rer[keep]
    # put phens in the same order as the rer values
    phens = phens[names(rer)]

    if((nb <- length(phens)) >= min.sp) {
      # for a binary or categorical trait, test that # species in fgd or per category is > min.pos
      if (method == "kw" || method == "aov") {
        counts = table(phens)
        if (length(counts) < 2 || min(counts) < min.pos) {
          next
        }
      }
      else if(method != "p" && sum(phens != 0) < min.pos) {
        next
      }

      x = rer
      if (!is.null(winsorizeRER)) {
        x = win(x, winsorizeRER)
      }
      if (!is.null(winsorizetrait)) {
        y = win(phens, winsorizetrait)
      }
      else {
        y = phens
      }
      if (method == "aov") {
        yfacts = as.factor(y)
        df = data.frame(x, yfacts)
        colnames(df) = c("RER", "category")
        ares = aov(RER ~ category, data = df)
        ares_Fval = summary(ares)[[1]][1, 4]
        ares_pval = summary(ares)[[1]][1, 5]
        corout[i, 1:3] = c(ares_Fval, nb, ares_pval)
        tukey = TukeyHSD(ares)
        groups = rownames(tukey[[1]])
        unnamedinds = which(is.na(names(tables)))
        if (length(unnamedinds > 0)) {
          newnamesinds = which(is.na(match(groups,
                                           names(tables))))
          if (length(newnamesinds) > 0) {
            names(tables)[unnamedinds][1:length(newnamesinds)] = groups[newnamesinds]
          }
        }
        for (k in 1:length(groups)) {
          name = groups[k]
          tables[[name]][i, "Rho"] = tukey[[1]][name,
                                                1]
          tables[[name]][i, "P"] = tukey[[1]][name,
                                              4]
        }
      }
      else if (method == "kw") {
        yfacts = factor(y)
        kres = kwdunn.test(x, yfacts, ncategories = lu)
        effect_size = kres$kw$H/(nb - 1)
        corout[i, 1:3] = c(effect_size, nb, kres$kw$p)
        
        for(k in 1:length(kres$dunn$Z)){ # length of kres$dunn$Z should be the same as length(tables) otherwise there's a problem
          tables[[k]][i, "Rho"] = kres$dunn$Z[k]
          tables[[k]][i, "P"] = kres$dunn$P.adjust[k]
        }
        
        # old code before speed up:
        # yfacts = as.factor(y)
        # df = data.frame(x, yfacts)
        # colnames(df) = c("RER", "category")
        # kres = kruskal.test(RER ~ category, data = df)
        # kres_Hval = kres$statistic
        # kres_pval = kres$p.value
        # corout[i, 1:3] = c(kres_Hval, length(phens), kres_pval)
        # dunn = dunnTest(RER ~ category, data = df,
        #                 method = "bonferroni")
        # groups = dunn$res$Comparison
        # unnamedinds = which(is.na(names(tables)))
        # if (length(unnamedinds > 0)) {
        #   newnamesinds = which(is.na(match(groups,
        #                                    names(tables))))
        #   if (length(newnamesinds) > 0) {
        #     names(tables)[unnamedinds][1:length(newnamesinds)] = groups[newnamesinds]
        #   }
        # }
        # for (k in 1:length(groups)) {
        #   name = groups[k]
        #   tables[[name]][i, "Rho"] = dunn$res$Z[k]
        #   tables[[name]][i, "P"] = dunn$res$P.adj[k]
        # }
      }
      else {
        cres = cor.test(x, y, method = method, exact = F)
        corout[i, 1:3] = c(cres$estimate, nb, cres$p.value)
      }
    }
  }

  # format and return the output
  corout = as.data.frame(corout)
  corout$p.adj = p.adjust(corout$P, method = "BH")
  if (method == "aov" || method == "kw") {
    for (i in 1:length(tables)) {
      # convert to a data frame
      tables[[i]] = as.data.frame(tables[[i]])
      # add an adjusted p value
      tables[[i]]$p.adj = p.adjust(tables[[i]]$P, method = "BH")
    }
    return(list(corout, tables))
  }
  else {
    corout
  }
}


#' main RER computation function
#' @param treesObj A treesObj created by \code{\link{readTrees}}
#' @param a cutoff value for branch lengths bellow which the branch lengths will be discarded, very data dependent but should roughly correspond to 0 or 1 sequence change on that branch. If left NULL this whill be set to the bottom 0.05 quantile. Set to 0 for no cutoff.
#' @param transform The transformation to apply to the trees branch values before computing relative rates. Available options are sqrt and log, sqrt is recommended.
#' @param weighted Use weighted regression to compute relative rates, meant to correct for the non-constant mean-variance relationship in evolutionary rate data.
#' @param useSpecies Give only a subset of the species to use for RER calculation. Some times excluding unusually long branches can provide more stable results
#' @param min.sp The minimum number of species needed to compute RER
#' @param scale Scale relative rates internally for each species subset. Increases computation time with little apparent benefit. Better to scale the final matrix.
#' @param doOnly The index of a specific tree in the treesObj to calculate RER for. Useful if a single result is needed quickly.
#' @param maxT The maximum number of trees to compute results for. Since this function takes some time this is useful for debugging.
#' @param plot Whether to plot the output of the correction for mean-variance relationship.
#' @return A numer of trees by number of paths matrix of relative evolutionary rates. Only an independent set of paths has non-NA values for each tree.
#' @export
getAllResiduals=function(treesObj, cutoff=NULL, transform="sqrt", weighted=T,  useSpecies=NULL,  min.sp=10, scale=T,  doOnly=NULL, maxT=NULL, scaleForPproj=F, mean.trim=0.05, plot=T){

  if(is.null(cutoff)){
    cutoff=quantile(treesObj$paths, 0.05, na.rm=T)
    message(paste("cutoff is set to", cutoff))
  }
  if (weighted){
    weights=computeWeightsAllVar(treesObj$paths, transform=transform, plot=plot)
    residfunc=fastLmResidMatWeighted
  }
  else{
    residfunc=fastLmResidMat
  }
  # residfunc=naresid

  if (is.null(useSpecies)){
    useSpecies=treesObj$masterTree$tip.label
    #mappedEdges=trees$mappedEdges
  }
  if(is.null(maxT)){
    maxT=treesObj$numTrees
  }
  if(transform!="none"){
    transform=match.arg(transform,c("sqrt", "log"))
    transform=get(transform)
  }
  else{
    transform=NULL
  }



  #cm is the names of species that are included in useSpecies and the master tree
  cm=intersect(treesObj$masterTree$tip.label, useSpecies)
  sp.miss = setdiff(treesObj$masterTree$tip.label, useSpecies)
  if (length(sp.miss) > 0) {
    message(paste0("Species from master tree not present in useSpecies: ", paste(sp.miss,
                                                                                 collapse = ",")))

  }

  rr=matrix(nrow=nrow(treesObj$paths), ncol=ncol(treesObj$paths))

  #maximum number of present species
  maxn=rowSums(treesObj$report[,cm])

  if(is.null(doOnly)){
    doOnly=1
  }
  else{
    maxT=1
  }
  skipped=double(nrow(rr))
  skipped[]=0

  for (i in doOnly:(doOnly+maxT-1)){

    if(sum(!is.na(rr[i,]))==0&&!skipped[i]==1){


      #get the ith tree
      tree1=treesObj$trees[[i]]

      #get the common species, prune and unroot
      both=intersect(tree1$tip.label, cm)
      if(length(both)<min.sp){
        next
      }
      tree1=unroot(pruneTree(tree1,both))

      #do the same for the refTree


      #find all the genes that contain all of the species in tree1
      allreport=treesObj$report[,both]
      ss=rowSums(allreport)
      iiboth=which(ss==length(both)) #this needs to be >1
      if (length(iiboth) < 2) {
        message(paste("Skipping i =",i,"(no other genes with same species set)"))
        next
      }

      nb=length(both)
      ai=which(maxn[iiboth]==nb)


      message(paste("i=", i))


      if(T){

        ee=edgeIndexRelativeMaster(tree1, treesObj$masterTree)

        ii= treesObj$matIndex[ee[, c(2,1)]]

        allbranch=treesObj$paths[iiboth,ii]
        if (is.null(dim(allbranch))) {
          message(paste("Issue with gettiing paths for genes with same species as tree",i))
          return(list("iiboth"=iiboth,"ii"=ii))
        }

        if(weighted){
          allbranchw=weights[iiboth,ii]
        }
        if(scaleForPproj){
          nv=apply(scaleMatMean(allbranch), 2, mean, na.rm=T, trim=mean.trim)
        }
        else{
          nv=apply(allbranch, 2, mean, na.rm=T, trim=mean.trim)
        }

        iibad=which(allbranch<cutoff)
        #don't scale
        #allbranch=scaleMat(allbranch)
        if(!is.null(transform)){
          nv=transform(nv)
          allbranch=transform(allbranch)
        }
        allbranch[iibad]=NA




        if(!scale){
          if(!weighted){
            proj=residfunc(allbranch[ai, ,drop=F], model.matrix(~1+nv))

          }
          else{

            proj=residfunc(allbranch[ai, ,drop=F], model.matrix(~1+nv), allbranchw[ai, ,drop=F])

          }
        }

        else{

          if(!weighted){
            proj=residfunc(allbranch[, ,drop=F], model.matrix(~1+nv))
          }
          else{

            proj=residfunc(allbranch[, ,drop=F], model.matrix(~1+nv),allbranchw)
          }

          proj=scale(proj, center = F)[ai, , drop=F]

        }


        #we have the projection



        rr[iiboth[ai],ii]=proj

      }

    }}
  message("Naming rows and columns of RER matrix")
  rownames(rr)=names(treesObj$trees)
  colnames(rr)=namePathsWSpecies(treesObj$masterTree)
  rr
}



#' turns a named vector of characters into a paths vector to be used with \code{\link{getAllCor}}
#' @param tip.vals the trait/phenotype/character value at the tip, \code{names(tip.vals)} should match some of the \code{mastertree$tip.label}, though a perfect match is not required
#' @param  treesObj A treesObj created by \code{\link{readTrees}}
#' @inheritParams  edgeVars
#' @return A vector of length equal to the number of paths in treesObj
#' @export
char2Paths=  function (tip.vals, treesObj, altMasterTree = NULL, metric = "diff",
                       se.filter = -1, ...)
{
  if (!is.null(altMasterTree)) {
    masterTree = altMasterTree
  }
  else if (!all(treesObj$masterTree$edge.length == 1)) {
    masterTree = treesObj$masterTree
  }
  else {
    message("The treesObj master tree has no edge lengths, please provide an alternative master tree")
    return()
  }
  cm=intersect(treesObj$masterTree$tip,intersect(names(tip.vals), masterTree$tip))

  #reduce to the same species set
  master.tree = pruneTree(masterTree, cm)
  tip.vals=tip.vals[cm]
  #make the tree with ancestral states
  charTree = edgeVars(master.tree, tip.vals, metric=metric, se.filter=se.filter, ...)


  sp.miss = setdiff(treesObj$masterTree$tip, names(tip.vals))
  if (length(sp.miss) > 0) {
    message(paste0("Species not present: ", paste(sp.miss,
                                                  collapse = ",")))

  }

  ap = allPaths(treesObj$masterTree)
  allPathMasterRelative(charTree, treesObj$masterTree, ap)
}





#' Creates paths from a set of foreground species
#' @param foreground. A character vector containing the foreground species
#' @param  treesObj A treesObj created by \code{\link{readTrees}}
#' @param plotTree Plot a tree representation of the result
#' @param clade A character string indicating which branches within the clade
#' containing the foreground species should be set to foreground. Must be one
#' of the strings "ancestral", "terminal", "all", or "weighted".
#' @param useSpecies Give only a subset of the species to use for ancestral state reconstruction
#' (e.g., only those species for which the trait can be reliably determined).
#' @return A vector of length equal to the number of paths in treesObj
#' @export
foreground2Paths = function(foreground,treesObj, plotTree=F, clade=c("ancestral","terminal","all","weighted"), useSpecies=NULL, transition="unidirectional"){
  #res = treesObj$masterTree
  #res$edge.length <- rep(0,length(res$edge.length))
  #res$edge.length[nameEdges(treesObj$masterTree) %in% foreground] = 1
  #names(res$edge.length) = nameEdges(treesObj$masterTree)
  res = foreground2Tree(foreground, treesObj, plotTree=plotTree, clade=clade, useSpecies=useSpecies, transition=transition)
  tree2Paths(res, treesObj)
}


#' Creates a binary trait tree from a set of foreground species.
#' @param foreground. A character vector containing the foreground species
#' @param treesObj A treesObj created by \code{\link{readTrees}}
#' @param collapse2anc Put all the weight on the ancestral branch when the trait appears on a whole clade
#' (redundant to "clade", kept for backwards compatibility)
#' @param plotTree Plot a tree representation of the result
#' @param wholeClade Whether to implement the weighted edge option across
#' all members of a foreground clade (redundant to "clade", kept for backwards compatibility)
#' @param clade A character string indicating which branches within the clade
#' containing the foreground species should be set to foreground. Must be one
#' of the strings "ancestral", "terminal", "all".
#' @param transition A character string indicating whether transitions between background and foreground branches
#' are "bidirectional" or "unidirectional" (no foreground to background transitions, the default)
#' @param weighted if set to TRUE weights foreground edges belonging to the same clade such that their branch lengths sum up to 1 (only done for clade options "all" and "terminal").
#' @param useSpecies Give only a subset of the species to use for ancestral state reconstruction
#' (e.g., only those species for which the trait can be reliably determined).
#' @return A tree with edge.lengths representing phenotypic states
#' @export
foreground2Tree = function(foreground,treesObj, plotTree=T, clade=c("ancestral","terminal","all"), weighted = F, transition = "unidirectional", useSpecies=NULL){
  clade <- match.arg(clade)
  res = treesObj$masterTree
  if (!is.null(useSpecies)) {
    sp.miss = setdiff(res$tip.label, useSpecies)
    if (length(sp.miss) > 0) {
      message(paste0("Species from master tree not present in useSpecies: ", paste(sp.miss,
                                                                                   collapse = ",")))
    }
    useSpecies = intersect(useSpecies, res$tip.label)
    res = pruneTree(res, useSpecies)
  } else {
    useSpecies = res$tip.label
  }
  foreground = intersect(foreground, useSpecies)
  res$edge.length <- rep(0,length(res$edge.length))
  if(clade == "terminal"){
    res$edge.length[nameEdges(res) %in% foreground] = 1
    names(res$edge.length) = nameEdges(res)
  }else if(clade == 'ancestral'){
    weighted = F
    if(transition == 'bidirectional'){
      res <- inferBidirectionalForegroundClades(res,foreground,ancestralOnly =T)
    }else{
      res <- inferUnidirectionalForegroundClades(res,foreground,ancestralOnly =T)
    }
  }else{
    if(transition == 'bidirectional'){
      res <- inferBidirectionalForegroundClades(res,foreground,ancestralOnly = F)
    }else{
      res <- inferUnidirectionalForegroundClades(res,foreground,ancestralOnly = F)
    }
  }
  if(weighted){
    if(clade == 'all'){
      tobeweighted <- rep(TRUE,length(res$edge.length))
      tobeweighted[res$edge.length == 0] <- FALSE
      while(sum(tobeweighted)>0){
        edgetodo <- which(tobeweighted == T)[1]
        clade.down.edges = getAllCladeEdges(res, edgetodo)
        if(length(clade.down.edges) > 1){
          clade.edges = c(clade.down.edges, edgetodo)
          clade.edges.toweight <- clade.edges[res$edge.length[clade.edges] == 1]
          res$edge.length[clade.edges.toweight] <- 1.0/(length(clade.edges.toweight))
          tobeweighted[clade.edges] <- FALSE
        } else{
          tobeweighted[clade.down.edges] <- FALSE
        }
      }
    } else if(clade == 'terminal'){
      tobeweightededgeterminalnode <- unique(res$edge[(res$edge[,2] %in% c(1:length(res$tip.label))),1])
      tobeweighted <- setdiff(match(tobeweightededgeterminalnode, res$edge[,2]), NA)
      for(edgetodo in tobeweighted){
        clade.down.edges = getAllCladeEdges(res, edgetodo)
        if(all(res$edge.length[clade.down.edges]==1)){
          res$edge.length[clade.down.edges] <- 0.5
        }
      }
    }
  }
  if(plotTree){
    res2=res
    mm=min(res2$edge.length[res2$edge.length>0])
    res2$edge.length[res2$edge.length==0]=max(0.02,mm/20)
    plot(res2, main = paste0("Clade: ",clade,'\nTransition: ',transition,'\nWeighted: ',weighted), cex = 0.5)
    if(weighted){
      labs <- round(res$edge.length,3)
      labs[labs == 0] <- NA
      edgelabels(labs, col = 'black', bg = 'transparent', adj = c(0.5,-0.5),cex = 0.4,frame='n')
    }
  }
  res
}


#'turns a named vector of characters into a paths vector to be used with \code{\link{getAllCor}} for categorical traits
#'@param tipvals the trait/phenotype/character value at the tip, \code{names(tip.vals)} should match some of the \code{mastertree$tip.label}, though a perfect match is not required
#'@param treesObj A treesObj created by \code{\link{readTrees}}
#'@param useSpecies Give only a subset of the species to use for ancestral state reconstruction
#' (e.g., only those species for which the trait can be reliably determined).
#'@param model Specifies what rate model to use
#'@param plot Plots a phenotype tree
#'@param anctrait The trait to use for all ancestral species instead of inferring ancestral states if not NULL. The default is NULL.
#'@param root_prior The prior probabilities of each trait at the root used to fit the transition matrix. Can be a vector of length equal to the number of states or one of the following: "flat", "empirical", "stationary", "likelihoods", "max_likelihood".
#'@return A vector of length equal to the number of paths in treesObj
#'@export
char2PathsCategorical = function(tipvals, treesObj, useSpecies = NULL,
                                 model = "ER", plot = FALSE, anctrait = NULL, root_prior = "auto") {
  #get tree
  tree = char2TreeCategorical(tipvals = tipvals, treesObj = treesObj,
                              useSpecies = useSpecies, model = model, plot = plot,
                              anctrait = anctrait, root_prior = root_prior)
  # get paths
  paths = tree2Paths(tree, treesObj, useSpecies = useSpecies,
                     categorical = T)
  return(paths)
}

#' the function root from ape modified to return additional information about node matching and which edges are inverted
#' ape citation: Paradis E. & Schliep K. 2019. ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R.Bioinformatics 35: 526-528.
#' @keywords internal
Root <- function(phy, outgroup, node = NULL, resolve.root = FALSE,
                 interactive = FALSE, edgelabel = FALSE, ...)
{
  if (!inherits(phy, "phylo"))
    stop('object not of class "phylo"')
  phy <- reorder(phy)
  n <- length(phy$tip.label)
  ROOT <- n + 1L

  if (interactive) {
    node <- identify(phy)$nodes
    cat("You have set resolve.root =", resolve.root, "\n")
  }

  ## added to solve some issues (2021-04-15):
  if (!interactive && is.null(node) && length(outgroup) > 1 && resolve.root)
    phy <- unroot(phy)
  ## -> the condition check should insure compatibility

  e1 <- phy$edge[, 1L]
  e2 <- phy$edge[, 2L]
  wbl <- !is.null(phy$edge.length)

  if (!is.null(node)) {
    if (node <= n) {
      print(n)
      stop("incorrect node#: should be greater than the number of taxa")
    }
    outgroup <- NULL
    newroot <- node
  } else {
    if (is.numeric(outgroup)) {
      if (any(outgroup > n))
        stop("incorrect taxa#: should not be greater than the number of taxa")
    }
    if (is.character(outgroup)) {
      outgroup <- match(outgroup, phy$tip.label)
      if (anyNA(outgroup))
        stop("specified outgroup not in labels of the tree")
    }
    if (length(outgroup) == n) return(phy)
    outgroup <- sort(outgroup) # used below

    ## First check that the outgroup is monophyletic, unless it has only one tip
    if (length(outgroup) > 1) {
      pp <- prop.part(phy)
      ingroup <- (1:n)[-outgroup]
      newroot <- 0L
      for (i in 2:phy$Nnode) {
        if (identical(pp[[i]], ingroup)) {
          ## inverted with the next if (... (2013-06-16)
          newroot <- e1[which(e2 == i + n)]
          break
        }
        if (identical(pp[[i]], outgroup)) {
          newroot <- i + n
          break
        }
      }
      if (!newroot)
        stop("the specified outgroup is not monophyletic")
      MRCA.outgroup <- i + n
    } else newroot <- e1[which(e2 == outgroup)]
  }

  N <- Nedge(phy)
  oldNnode <- phy$Nnode

  Nclade <- tabulate(e1)[ROOT] # degree of the root node
  ## if only 2 edges connect to the root, we have to fuse them:
  fuseRoot <- Nclade == 2

  if (newroot == ROOT) {
    if (!resolve.root) return(phy) # else (resolve.root == TRUE)
    if (length(outgroup) > 1) outgroup <- MRCA.outgroup
    if (!is.null(node))
      stop("ambiguous resolution of the root node: please specify an explicit outgroup")

    k <- which(e1 == ROOT) # find the basal edges
    if (length(k) > 2) {
      i <- which(e2 == outgroup) # outgroup is always of length 1 here
      j <- k[k != i]
      newnod <- oldNnode + n + 1L
      phy$edge[j, 1] <- newnod

      phy$edge <- rbind(c(ROOT, newnod), phy$edge)
      if (wbl) phy$edge.length <- c(0, phy$edge.length)

      phy$Nnode <- phy$Nnode + 1L
    }
  } else {
    phy$root.edge <- NULL # just in case

    INV <- logical(N)
    w <- which(e2 == newroot)
    anc <- e1[w]
    i <- w

    nod <- anc

    if (nod != ROOT) {
      INV[w] <- TRUE
      i <- w - 1L
      repeat {
        if (e2[i] == nod) {
          if (e1[i] == ROOT) break
          INV[i] <- TRUE
          nod <- e1[i]
        }
        i <- i - 1L
      }
    }

    ## we keep the edge leading to the old root if needed:
    if (!fuseRoot) INV[i] <- TRUE

    ## bind the other clades...

    if (fuseRoot) { # do we have to fuse the two basal edges?
      k <- which(e1 == ROOT)
      k <- if (k[2] > w) k[2] else k[1]
      phy$edge[k, 1] <- phy$edge[i, 2]
      if (wbl)
        phy$edge.length[k] <- phy$edge.length[k] + phy$edge.length[i]
    }

    if (fuseRoot) phy$Nnode <- oldNnode - 1L

    ## added after discussion with Jaime Huerta Cepas (2016-07-30):
    if (edgelabel) {
      phy$node.label[e1[INV] - n] <- phy$node.label[e2[INV] - n]
      phy$node.label[newroot - n] <- ""
    }

    phy$edge[INV, ] <- phy$edge[INV, 2:1]

    if (fuseRoot) {
      phy$edge <- phy$edge[-i, ]
      if (wbl) phy$edge.length <- phy$edge.length[-i]
      N <- N - 1L
    }

    if (resolve.root) {
      newnod <- oldNnode + n + 1L
      if (length(outgroup) == 1L) {
        wh <- which(phy$edge[, 2] == outgroup)
        #phy$edge[1] <- newnod
        k <- which(phy$edge[, 1] == newroot) # wh should be among k
        phy$edge[k[k != wh], 1] <- newnod
        o <- c((1:N)[-wh], wh)
        phy$edge <- rbind(c(newroot, newnod), phy$edge[o, ])
        if (wbl) phy$edge.length <- c(0, phy$edge.length[o])
      } else {
        wh <- which(phy$edge[, 1] == newroot)
        phy$edge[wh[-1], 1] <- newnod
        s1 <- 1:(wh[2] - 1)
        s2 <- wh[2]:N
        phy$edge <-
          rbind(phy$edge[s1, ], c(newroot, newnod), phy$edge[s2, ])
        if (wbl)
          phy$edge.length <- c(phy$edge.length[s1], 0, phy$edge.length[s2])
      }
      phy$Nnode <- phy$Nnode + 1L
    }
  }
  ## The block below renumbers the nodes so that they conform
  ## to the "phylo" format
  newNb <- integer(n + phy$Nnode)
  newNb[newroot] <- n + 1L
  sndcol <- phy$edge[, 2] > n
  newNb[sort(phy$edge[sndcol, 2])] <- n + 2:phy$Nnode
  phy$edge[sndcol, 2] <- newNb[phy$edge[sndcol, 2]]
  phy$edge[, 1] <- newNb[phy$edge[, 1]]

  if (!is.null(phy$node.label)) {
    newNb <- newNb[-(1:n)]
    if (fuseRoot) {
      newNb <- newNb[-1]
      phy$node.label <- phy$node.label[-1]
    }
    phy$node.label <- phy$node.label[order(newNb)]
    if (resolve.root) {
      phy$node.label[is.na(phy$node.label)] <- phy$node.label[1]
      phy$node.label[1] <- "Root"
    }
  }
  attr(phy, "order") <- NULL
  reorder.phylo(phy)

  # make a map of old node numbers to new node numbers
  map = cbind((n + 1L):(n + oldNnode), newNb[(n + 1L):(n + oldNnode)])
  colnames(map) = c("old", "new")
  # return the extra info
  return(list(phy = phy, map = map, flipped_edges = INV))
}

#' Calculates ancestral likelihoods of ancestral states for discrete characters
#' @param tree An object of class "phylo"
#' @param tips The phenotype data for the species at the tips of the tree. Must be in the same order as tree$tip.label
#' @param Q A transition rate matrix fit on the tree and tips
#' @return A matrix of ancestral likelihoods with rows in the same order as the internal nodes in the tree
#' @export
asymmRerootingMethod <- function(tree, tips, Q) {

  ntips = length(tree$tip.label)
  # matrix to store the likelihoods
  marg_anc_liks = matrix(nrow = tree$Nnode, ncol = nrow(Q))
  # convert tips to a matrix
  tips = to.matrix(tips, sort(unique(tips)))

  # loop through the internal nodes
  nodes = unique(tree$edge[,1]) # get the internal nodes

  for(n in nodes) {
    # reroot at node n
    trInfo = Root(tree, node = n)

    # check if it's the original tree
    if(inherits(trInfo, "phylo")) {
      # this is the original tree
      tr = reorder(trInfo, order = "postorder")

      # matrix to store likelihoods during pruning algorithm
      liks = matrix(nrow = tr$Nnode, ncol = nrow(Q))
      # add tips to liks
      liks = rbind(tips, liks)

      # loop through the internal nodes
      parents = unique(tr$edge[,1])
      for(i in 1:length(parents)) {
        # get children of parent node
        p = parents[i]
        cc = tr$edge[,2][which(tr$edge[,1] == p)] # children nodes
        ee = tr$edge.length[which(tr$edge[,1] == p)] # edge lengths

        v = vector(mode = "list", length = length(cc))
        for(c in 1:length(cc)) {
          P = expm(Q * ee[c])
          v[[c]] = P %*% liks[cc[c],]
        }
        ll = Reduce("*",v)[,1]
        liks[p,] =  ll/sum(ll)  # normalize and store in liks
      }

      # add the root likelihood to marg_anc_liks
      marg_anc_liks[(n - ntips),] = liks[n,] # n is the root

    }
    else { # not the original tree
      # get the re-rooted tree
      tr = trInfo$phy
      map = trInfo$map

      # reorder - postorder traversal for pruning algorihtm
      tr = reorder(tr, order = "postorder")

      # matrix to store likelihoods during pruning algorithm
      liks = matrix(nrow = tr$Nnode, ncol = nrow(Q))
      # add tips to liks
      liks = rbind(tips, liks)

      # loop through the internal nodes
      parents = unique(tr$edge[,1])
      for(i in 1:length(parents)) {
        # get children of parent node
        p = parents[i]
        oldP = map[,1][which(map[,2] == p)]

        cc = tr$edge[,2][which(tr$edge[,1] == p)] # children nodes
        oldcc = sapply(cc, function(x){
          if(x <= ntips) x
          else map[,1][which(map[,2] == x)]
        })

        ee = tr$edge.length[which(tr$edge[,1] == p)] # edge lengths

        v = vector(mode = "list", length = length(cc))
        for(c in 1:length(cc)) {
          # find the edge defined by oldP and oldcc[c]
          oldEdge = which(apply(tree$edge,1,function(x){oldP %in% x && oldcc[c] %in% x}))

          # determine if flipped or not flipped
          if(length(oldEdge > 0)) { # if the oldEdge exists
            flipped = trInfo$flipped_edges[oldEdge]
          } else {
            flipped = FALSE
          }

          # if flipped reverse the dot product between P and liks
          P = expm(Q * ee[c])
          if(flipped) {
            v[[c]] = t(liks[cc[c],] %*% P) # take transpose so it has one column with multiple rows
          } else {
            v[[c]] = P %*% liks[cc[c],]
          }
        }
        ll = Reduce("*",v)[,1]
        liks[p,] =  ll/sum(ll)  # normalize and store in liks
      }

      # add the root likelihood to marg_anc_liks
      root = min(parents)
      marg_anc_liks[(n - ntips),] = liks[root,] # n is the root in the old tree
    }
  }
  return(marg_anc_liks)
}

#' Returns ancestral likelihoods at each node. Based on code from ace in ape and fitMk/rerootingMethod in phytools
#' @param tree Object of class phylo that has been pruned to only contain the species in tipvals
#' @param tipvals The phenotype data for the extant species in the tree in the same order as tree$tip.label and mapped to integers
#' @param Q A transition matrix, if NULL the transition matrix is fit with fit_mk from castor package
#' @param rate_model The rate model to use for fitting the transition matrix if one is not provided
#' @param root_prior The root prior used when fitting the transition matrix if one is not provided
#' @return The ancestral likelihoods at each node in the tree
#' @export
getAncLiks <- function(tree, tipvals, Q = NULL, rate_model = "ER", root_prior = "auto") {

  ntips = length(tree$tip.label) # number of tips
  nstates = length(unique(tipvals)) # number of states (categories)

  # matrix to store the likelihoods
  liks = matrix(nrow = tree$Nnode, ncol = nstates)
  # convert tips to a matrix
  tips = to.matrix(tipvals, sort(unique(tipvals)))
  liks = rbind(tips, liks)

  # get transition matrix, Q
  if(is.null(Q)) {
    # intlabels = map_to_state_space(tipvals)
    # Q = fit_mk(trees = tree, Nstates = intlabels$Nstates,
    #            tip_states = intlabels$mapped_states,
    #            rate_model = rate_model, root_prior = root_prior)$transition_matrix
    Q = fit_mk(trees = tree, Nstates = length(unique(tipvals)),
               tip_states = tipvals,
               rate_model = rate_model, root_prior = root_prior)$transition_matrix
  }

  tree = reorder(tree, order = "postorder")

  # forward pass:
  parents = unique(tree$edge[,1])
  for(i in 1:length(parents)) {
    # get children of parent node
    p = parents[i]
    cc = tree$edge[,2][which(tree$edge[,1] == p)] # children nodes
    ee = tree$edge.length[which(tree$edge[,1] == p)] # edge lengths

    v = vector(mode = "list", length = length(cc))
    for(c in 1:length(cc)) {
      P = expm(Q * ee[c])
      v[[c]] = P %*% liks[cc[c],]
    }
    ll = Reduce("*",v)[,1]
    liks[p,] =  ll/sum(ll)  # normalize and store in liks
  }

  # backward pass:
  for(i in length(parents):1){
    # get children of parent node
    p = parents[i]
    cc = tree$edge[,2][which(tree$edge[,1] == p)] # children nodes
    ee = tree$edge.length[which(tree$edge[,1] == p)] # edge lengths

    for(c in 1:length(cc)){
      des = cc[c]
      if(des > ntips) { # if the child node is an internal node, update likelihood
        P = expm(Q * ee[c])
        tmp <- Matrix::t(liks[p, ] / (P %*% liks[des, ])) # from ace, but flip dot product and take transpose
        ll <- (tmp %*% P) * liks[des, ] # from ace
        ll = ll[1,]
        liks[des,] = ll/sum(ll, na.rm = TRUE)
      }
    }
  }
  liks = liks[-(1:ntips),]
  if(sum(is.nan(liks)) > 0) {
    warning("NaN values produced. Zero values in the transition rate matrix preventing necessary transitions. Use a different rate model.")
  }
  return(liks)
}


#' Creates a categorical trait tree from a set of tip species.
#'@param tipvals the trait/phenotype/character value at the tip, \code{names(tip.vals)} should match some of the \code{mastertree$tip.label}, though a perfect match is not required
#'@param treesObj A treesObj created by \code{\link{readTrees}}
#'@param useSpecies Give only a subset of the species to use for ancestral state reconstruction
#' (e.g., only those species for which the trait can be reliably determined).
#'@param model Specifies what rate model to use
#'@param plot Plots a phenotype tree
#'@param anctrait The trait to use for all ancestral species instead of inferring ancestral states if not NULL. The default is NULL.
#'@param root_prior The prior probabilities of each trait at the root used to fit the transition matrix. Can be a vector of length equal to the number of states or one of the following: "flat", "empirical", "stationary", "likelihoods", "max_likelihood".
#'@return A tree with edge.length representing phenotype states
#'@export
char2TreeCategorical <- function (tipvals, treesObj, useSpecies = NULL,
                                     model = "ER", root_prior = "auto",
                                     plot = FALSE, anctrait = NULL)
{
  # get the master tree and prune to include useSpecies/species with phenotype data
  mastertree = treesObj$masterTree
  if (!is.null(useSpecies)) {
    sp.miss = setdiff(mastertree$tip.label, useSpecies)
    if (length(sp.miss) > 0) {
      message(paste0("Species from master tree not present in useSpecies: ",
                     paste(sp.miss, collapse = ",")))
    }
    useSpecies = intersect(mastertree$tip.label, useSpecies)
    mastertree = pruneTree(mastertree, useSpecies)
    # unroot the tree after pruning
    mastertree = unroot(mastertree)
  }
  else {
    mastertree = pruneTree(mastertree, intersect(mastertree$tip.label,
                                                 names(tipvals)))
    # unroot the tree after pruning
    mastertree = unroot(mastertree)
  }
  # use ASR to infer phenotype tree
  if (is.null(anctrait)) {

    tipvals <- tipvals[mastertree$tip.label]
    intlabels <- map_to_state_space(tipvals)
    print("The integer labels corresponding to each category are:")
    print(intlabels$name2index)

    ancliks = getAncLiks(mastertree, intlabels$mapped_states, rate_model = model,
                         root_prior = root_prior)

    states = rep(0, nrow(ancliks))
    for (i in 1:length(states)) {
      states[i] = which.max(ancliks[i,])
    }
    states = c(intlabels$mapped_states, states)
    tree = mastertree
    tree$edge.length = states[tree$edge[, 2]]

    # convert to binary tree if necessary, plot, & return tree
    if(length(unique(tipvals)) == 2) {
      if(sum(! unique(tipvals) %in% c(TRUE,FALSE)) > 0) { # check that the two categories are TRUE/FALSE
        message("Returning categorical tree for binary phenotype because phenotype values are not TRUE/FALSE")
      } else {
        tree$edge.length = ifelse(tree$edge.length == 2, 1, 0)
        print("There are only 2 categories: returning a binary phenotype tree.")
        if (plot) {
          plotTree(tree)
        }
        return(tree)
      }
    }
    if (plot) {
      plotTreeCategorical(tree, category_names = intlabels$state_names,
                          master = mastertree, node_states = states)
    }
    return(tree)
  }
  else {
    if (length(unique(tipvals)) <= 2) {
      fgspecs <- names(tipvals)[tipvals != anctrait]
      res <- foreground2Tree(fgspecs, treesObj, plotTree = plot,
                             clade = "terminal", useSpecies = useSpecies)
      print("There are only 2 categories: returning a binary phenotype tree.")
      if(plot) {
        plotTree(res)
      }
      return(res)
    }
    else {
      tipvals <- tipvals[mastertree$tip.label]
      intlabels <- map_to_state_space(tipvals)
      j <- which(intlabels$state_names == anctrait)
      if (length(j) < 1) {
        warning("The ancestral trait provided must match one of the traits in the phenotype vector.")
      }
      res = mastertree
      res$edge.length <- rep(j, length(res$edge.length))
      traits <- intlabels$state_names
      for (trait in traits) {
        if (trait == anctrait) {
          next
        }
        i <- which(intlabels$state_names == trait)
        res$edge.length[nameEdges(res) %in% names(tipvals)[tipvals ==
                                                             trait]] = i
      }
      names(res$edge.length) = nameEdges(res)
      if (plot) {
        # get states for plotting
        states = res$edge.length[order(res$edge[,2])]
        states = c(j, states) # add root since it's not included in res$edge.length (no edge leading to the root)
        plotTreeCategorical(res, category_names = traits,
                            master = treesObj$masterTree,
                            node_states = states)
      }
      print("Category names are mapped to integers as follows:")
      print(intlabels$name2index)
      return(res)
    }
  }
}

#' Returns a list of states in node order from a phenotype tree
#' @param phenTree The phenotype tree returned by char2TreeCategorical
#' @param root_state The state of the root node, this is not stored along the edge lengths of the tree
#' @param tipvals The states of the tips of the tree in order of phenTree$tip.label
#' @return Returns a vector of the category states of species in order of the nodes in the tree
#' @export
getStatesFromPhenTree <- function(phenTree, root_state, tipvals) {
  ntips = length(phenTree$tip.label)
  node_states = phenTree$edge.length[order(phenTree$edge[,2])][(ntips + 1):(phenTree$Nnode + ntips - 1)]
  states = c(tipvals, root_state, node_states)
  return(states)
}

#' @keywords internal
inferUnidirectionalForegroundClades <- function(tree, fgd = NULL, ancestralOnly = F){
  finaltree <- tree
  finaltree$edge.length <- rep(0, length(tree$edge.length))
  finaltree$edge.length[nameEdges(finaltree) %in% fgd] <- 1
  #figure out node depth - terminal nodes have depth of 1; higher numbers indicate ancestral nodes;
  nodedepths <- node.depth(finaltree)
  edgeterminalnodedepths <- nodedepths[finaltree$edge[,2]]
  #going from 1-away from terminal ancestral branch to the base of the tree, figure out branches where all downstream clades are foreground
  for(inode in sort(unique(edgeterminalnodedepths))[-1]){
    edgesToDo <- which(edgeterminalnodedepths == inode)
    for(edgeindex in edgesToDo){
      clade.edges = getAllCladeEdges(finaltree, edgeindex)
      if(all(finaltree$edge.length[clade.edges]==1)){
        finaltree$edge.length[edgeindex] <- 1
      }
    }
  }
  if(ancestralOnly){
    for(edgeii in 1:length(finaltree$edge.length)){
      if(finaltree$edge.length[edgeii] == 1){
        if(nameEdges(finaltree)[edgeii]==""){
          clade.edges = setdiff(getAllCladeEdges(finaltree, edgeii), edgeii)
          finaltree$edge.length[clade.edges] <- 0
        }
      }
    }
  }
  finaltree
}

#' @keywords internal
inferBidirectionalForegroundClades <- function(treeinput, foreground = NULL, ancestralOnly = F){
  tree <- treeinput
  tip.vals=rep(0, length(tree$tip.label))
  names(tip.vals)=tree$tip.label
  tip.vals[foreground]=1
  tmp=cbind(as.character(tip.vals))
  rownames(tmp)=names(tip.vals)
  tip.vals=tmp
  #Add option to function for "type" within ancestral.pars
  ancres=ancestral.pars(tree, df<-as.phyDat(tip.vals, type="USER", levels=unique(as.character(tip.vals))),type="ACCTRAN" )
  ancres=unlist(lapply(ancres, function(x){x[2]}))
  internalVals=ancres
  #evals=matrix(nrow=nrow(treesObj$masterTree$edge), ncol=2)
  evals=matrix(nrow=nrow(tree$edge), ncol=2)
  eres=ancres
  #evals[,1]=eres[treesObj$masterTree$edge[,1]]
  evals[,1]=eres[tree$edge[,1]]
  #evals[,2]=eres[treesObj$masterTree$edge[,2]]
  evals[,2]=eres[tree$edge[,2]]
  tree$edge.length=evals[,2]-evals[,1]
  #res$edge.length[res$edge.length<1]=0
  if(!ancestralOnly){
    edgeIndex=which(tree$edge.length>0)
    edgeIndexNeg=which(tree$edge.length<0)
    edgeIndexAll = c(edgeIndex,edgeIndexNeg)
    edgeDirection = c(rep(1, length(edgeIndex)),rep(-1, length(edgeIndexNeg)))
    edgedf = data.frame(edgeIndexAll,edgeDirection)
    edgedf = edgedf[order(edgedf$edgeIndexAll),]
    clade.edges=NA
    clade.lengths=NA
    cladedf = data.frame(clade.edges,clade.lengths)
    for(i in 1:nrow(edgedf)) { #Does this go from ancestral to terminal?
      #save the clade until the edges no longer overlap
      clade.edges=getAllCladeEdges(tree, edgedf$edgeIndexAll[i])
      clade.edges=unique(c(edgedf$edgeIndexAll[i], clade.edges))
      if (any(clade.edges %in% cladedf$clade.edges)==F) {
        tree$edge.length[cladedf$clade.edges[which(cladedf$clade.lengths==1)]]=1
        if (edgedf$edgeDirection[i] == 1) {
          clade.lengths = c(rep(1,length(clade.edges)))
        } else {
          clade.lengths = c(rep(0,length(clade.edges)))
        }
        cladedf = data.frame(clade.edges,clade.lengths)
      } else {
        #update df lengths
        if (edgedf$edgeDirection[i] == 1) {
          cladedf$clade.lengths[which(cladedf$clade.edges %in% clade.edges)] = 1
        } else {
          cladedf$clade.lengths[which(cladedf$clade.edges %in% clade.edges)] = 0
        }
      }
    }
    #update edge lengths from the final clade
    tree$edge.length[cladedf$clade.edges[which(cladedf$clade.lengths==1)]]=1
    tree$edge.length[tree$edge.length<0]=0
  }
  tree$edge.length[tree$edge.length<0]=0
  tree
}

#' @keywords internal
nameEdges=function(tree){
  nn=character(nrow(tree$edge))
  iim=match(1:length(tree$tip.label), tree$edge[,2])
  nn[iim]=tree$tip.label
  nn
}

#' pseudoroot trait tree to match the psuedoroot of trees from \code{\link{readTrees}}  if the trees are reconcilable
#' @param tree A trait tree with branch lengths representing trait values
#' @param treesObj A treesObject created by \code{\link{readTrees}}
#' @return A trait tree with the correct topology
#' @export
fixPseudoroot=function(tree, treesObj){
  if(RF.dist(tree, treesObj$masterTree)>0){
    stop("The trait tree and treesObj$masterTree are not reconcilable - they have different topologies")
  }
  #fix pseudorooting
  tr1=tree
  tr2=treesObj$masterTree
  #get species at pseudoroot
  toroot=tr2$tip.label[tr2$edge[,2][tr2$edge[,1]==as.numeric(names(which(table(tr2$edge[,1])==3)))]]
  toroot=toroot[!is.na(toroot)]
  #pick one, must be in tr1
  if(toroot[[1]] %in% tr1$tip.label){
    toroot=toroot[[1]]
  }else if(toroot [[2]] %in% tr1$tip.label){
    toroot=toroot[[2]]
  }else{
    stop("Key species missing from trait tree")
  }
  tr1=root.phylo(tr1, toroot)
  tree=tr1
  plot(tree)
  message("Make sure the branch lengths for the new trait tree are correct")
  return(tree)
}

#' Generate a phenotype paths vector from a phenotype tree
#'
#' \code{tree2Paths} generates a phenotype paths vector matching the treesObject
#'     from a tree where branches specify phenotypes.
#'
#' The tree topology of the phenotype tree must match that of the master tree within the treesObject.
#'
#' @param tree A phenotype tree, with branch length encoding a phenotype.
#' @param treesObj A treesObject created by \code{\link{readTrees}}
#' @param binarize Force binary path representation. Default action depends upon the type of data within the phenotype tree
#'     (binary or continuous).
#'     \itemize{
#'    \item If binary (all branch lengths == 0 or 1): Sets all positive path values to 1. Useful if the tree has non-zero branch lengths
#'        for an internal branch or branches; otherwise, values are simply added along branches when calculating paths.
#'        Default behavior: binarize = TRUE.
#'    \item If continuous (not all branch lengths == 0 or 1): Sets all path values > the mean to 1 and all those <= the mean to 0.
#'        Converts a continuous phenotype to a binary phenotype, with state determined by comparison to the mean across all paths.
#'        Default behavior: binarize = FALSE.
#'        }
#' @param useSpecies Give only a subset of the species to use for ancestral state reconstruction
#' (e.g., only those species for which the trait can be reliably determined).
#' @return A vector of length equal to the number of paths in treesObj
#' @export
tree2Paths=function(tree, treesObj, binarize=NULL, useSpecies=NULL, categorical = F){
  stopifnot(class(tree)[1]=="phylo")
  stopifnot(class(treesObj)[2]=="treesObj")

  isbinarypheno <- sum(tree$edge.length %in% c(0,1)) == length(tree$edge.length) #Is the phenotype tree binary or continuous?
  if (is.null(binarize)) { #unless specified, determine default for binarize based on type of phenotype tree
    if (isbinarypheno) {
      binarize = T #default for binary phenotype trees: set all positive paths = 1
    } else {
      binarize = F #default for continuous phenotype trees: do not convert to binary
    }
  }
  #unroot if rooted
  if (is.rooted(tree)) {
    tree = unroot(tree)
  }

  #reduce tree to species in master tree and useSpecies
  sp.miss = setdiff(tree$tip.label, union(treesObj$masterTree$tip.label, useSpecies))
  if (length(sp.miss) > 0) {
    message(paste0("Species from tree not present in master tree or useSpecies: ", paste(sp.miss,
                                                                                         collapse = ",")))
  }
  if (!is.null(useSpecies)) {
    tree = pruneTree(tree, intersect(intersect(tree$tip.label, treesObj$masterTree$tip.label), useSpecies))
  } else {
    tree = pruneTree(tree, intersect(tree$tip.label, treesObj$masterTree$tip.label))
  }

  treePaths=allPaths(tree, categorical = categorical)
  map=matchAllNodes(tree,treesObj$masterTree)

  #remap the nodes
  treePaths$nodeId[,1]=map[treePaths$nodeId[,1],2 ]
  treePaths$nodeId[,2]=map[treePaths$nodeId[,2],2 ]

  #indices for which paths to return
  ii=treesObj$ap$matIndex[(treePaths$nodeId[,2]-1)*nrow(treesObj$ap$matIndex)+treePaths$nodeId[,1]]

  vals=double(length(treesObj$ap$dist))
  vals[]=NA
  vals[ii]=treePaths$dist
  if(binarize){
    if(isbinarypheno) {
      vals[vals>0]=1
    } else {
      mm=mean(vals)
      vals[vals>mm]=1
      vals[vals<=mm]=0
    }
  }
  vals
}

#' Makes a binary path vector from either a tree of class "phylo" or a foreground species set supplied as a character vector
#' @param input Either a phenotype tree of class "phylo" (with branch length encoding a phenotype) or a character vector of foreground branches
#' @param  treesObj A treesObj created by \code{\link{readTrees}}
#' @return A vector of length equal to the number of paths in treesObj
#' @export
makeBinaryPaths=function(input, treesObj){
  if(class(input)=="character"){
    foreground2Paths(input, treesObj)
  }
  else if (class(input)=="phylo"){
    tree2Paths(input, treesObj)
  }
  else{
    message("Need either a phylo tree or a character vector")
    return(NULL)
  }
}


#' maps a vector of traits onto a reference tree
#' @param mastertree the tree species the topology of the output tree and the branch lengths are used to infer ancestral states
#' @param tip.vals the trait/phenotype/character value at the tip, \code{names(tip.vals)} should match some of the \code{mastertree$tip.label}, though a perfect match is not required
#' @param metric The metric used to translate node values into branch values. "Diff" takes the difference and makes the result phylogenetically independent. Other possible values are "mean" (the mean of the two values) and "last" the value of the most recent species on the branch. The last two options are not phylogenetically independent and downstream computations for those are not yet implemented
#' @param se.filter Will remove branch values that are not at least \code{se.filter*edge.se} away from 0 (where edge.se is the standard error in the estimate for the edge value). Only implemented for \code{metric="diff"}. Set \code{se.filter} to a positive value to filter. By default no filtering is done.
#' @param return.var Returns the variance instead of the mean. Useful for seeing which estimates have high confidence.
#' @return A phylo tree with branch values computed from the input tip.values
#' @export
edgeVars=function(mastertree,tip.vals, metric="diff", se.filter=-1, return.var=F){
  message(paste0("using metric ", metric, ", with filtering constant ", se.filter))
  metric=match.arg(metric, c("diff", "mean", "last"))
  cm=intersect(mastertree$tip.label, names(tip.vals))
  mastertree=pruneTree(mastertree, cm)

  #sets edge length to the difference between two nodes, the evolutionary change (i.e. the change between species A and its ancestral species)
  tip.vals=tip.vals[mastertree$tip.label]
  res=fastAnc(mastertree, x=tip.vals, vars=T)
  vars=c(rep(NA, length(tip.vals)), res$var)
  res=c(tip.vals, res$ace)
  names(res)[1:length(tip.vals)]=as.character(1:length(tip.vals))
  evals=matrix(nrow=nrow(mastertree$edge), ncol=2)
  evars=matrix(nrow=nrow(mastertree$edge), ncol=2)

  evals[,1]=res[mastertree$edge[,1]]
  evals[,2]=res[mastertree$edge[,2]]
  evars[,1]=vars[mastertree$edge[,1]]
  evars[,2]=vars[mastertree$edge[,2]]
  newtree=mastertree

  if(metric=="diff"){
    newtree$edge.length=evals[,2]-evals[,1]
    edge.se=sqrt(apply(evars,1,mean, na.rm=T))/sqrt(length(cm))
    #   plot(newtree$edge.length, edge.se)
    iibad=which(abs(newtree$edge.length)<se.filter*edge.se)
    newtree$edge.length[iibad]=NA
  }
  else if (metric=="mean"){
    newtree$edge.length=evals[,2]+evals[,1]
  }
  else if (metric=="last"){
    newtree$edge.length=evals[,1]
    evars[,2]=0
  }
  if(!return.var){
    return(newtree)
  }
  else{
    message("Returning variance")
    newtree$edge.length=(edge.se*length(cm))^2

    return(newtree)
  }
}


#' wrapper around \code{\link[ape]{drop.tip}}
#' @param  tree a "phylo" tree
#' @param  tip.names The tip names to keep in the tree
#' @return  A new pruned tree
#' @export
pruneTree=function(tree, tip.names){
  keep=intersect(tree$tip.label, tip.names)
  torm=setdiff(tree$tip.label, keep)
  tree=drop.tip(tree, torm)
  tree
}



#' @keywords  internal
transformMat=function(tree){

  nA=length(tree$tip.label)+tree$Nnode
  matIndex=matrix(nrow=nA, ncol=nA)
  mat=matrix(nrow=nrow(tree$edge), ncol=0)
  index=1
  for ( i in 1:nA){
    ia=getAncestors(tree,i)
    if(length(ia)>0){
      thisindex=double()
      ia=c(i,ia)
      for (k in 2:length(ia)){
        j=ia[k]

        thisindex=c(thisindex, which(tree$edge[,2]==ia[k-1]&tree$edge[,1]==ia[k]))

        vals=rep(0, nrow(mat))
        vals[thisindex]=1
        mat=cbind(mat, vals)
      }
    }
  }
  mat
}







edgeIndexRelativeMaster=function(tree, masterTree){
  map=matchAllNodes(tree,masterTree)
  newedge=tree$edge
  newedge[,1]=map[newedge[,1],2]
  newedge[,2]=map[newedge[,2],2]
  newedge
}

edgeReOrder=function(tree, masterTree){
  map=matchAllNodes(tree,masterTree)
  #rename the tree edge
  newedge=tree$edge
  newedge[,1]=map[newedge[,1],2]
  newedge[,2]=map[newedge[,2],2]
  edgenames=namePaths(newedge)

  masternames=namePaths(masterTree$edge)

  match(masternames, edgenames)
}

namePaths=function(nodeMat, invert=F, mult=1000){
  warning("Why am I doing this")
  if(invert){
    nodeMat=nodeMat[,c(2,1)]
  }
  return(nodeMat[,1]*mult+nodeMat[,2])
}
printTipDist=function(tree, node){
  if(is.character(node)){
    node=match(node, tree$tip.label)
  }
  ii=which(tree$edge[,2]==node)
  print(tree$edge.length[ii])
}
printLeaveEdges=function(tree){
  nA=tree$Nnode+length(tree$tip)
  nT=length(tree$tip.label)
  ii=match(1:nT, tree$edge[,2])
  tmp=as.data.frame(tree$edge)
  tmp[ii,2]=tree$tip.label[tmp[ii,2]]
  show(tmp[ii,])
  tmp
}
printClade=function(tree,node){
  nT=length(tree$tip.label)
  clade=character()
  if(node<=nT){
    clade=c(clade, tree$tip.label[node])
  }
  else{
    ic=getChildren(tree,node)
    for ( i in 1:length(ic)){
      clade=c(clade, printClade(tree,ic[i]))
    }
  }
  paste(clade, collapse=",")
}

CanonicalForm=function(tree){
  par(mfrow=c(1,2))
  #  plot(tree)
  oo=order(tree$tip.label)
  tree$tip.label=tree$tip.label[oo]
  ii=match(1:length(oo), tree$edge[,2])
  tree$edge[ii,2]=order(oo)
  #plot(tree)
  rotateConstr(tree, sort(tree$tip.label))

}



rescaleTree=function(tree){
  tree$edge.length=tree$edge.length/sqrt(sum(tree$edge.length^2))
  tree
}
treeSum=function(tree){
  sum(tree$edge.length^2)
}
distToVec=function(dist){
  vec=as.vector(dist)
  names=character()
  for(i in 1:nrow(dist)){
    for(j in 1:ncol(dist)){
      names=c(names, paste(rownames(dist)[i], colnames(dist)[j], sep="_"))
    }
  }
  names(vec)=names
  vec
}






#linear fit with intercept
residLN=function(x,y, plot=F){
  if(plot){
    plot(y,x);abline(a=0,b=1)
  }
  return(as.vector(resid(x, cbind(rep(1,length(y)), y))))
}
#linear fit no intersept
residLN0=function(x,y, plot=F){
  if(plot){
    plot(y,x);abline(a=0,b=1)
  }
  return(as.vector(resid(x, cbind(y))))
}


#linear fit in sqroot space
residSQ=function(x,y, plot=F){
  x=sqrt(x);y=sqrt(y)
  if(plot){
    plot(y,x);abline(a=0,b=1)
  }
  return(as.vector(resid(x, cbind(rep(1,length(y)), y))))
}


#loess fit
residLO=function(x,y, plot=F){
  x=as.vector(sqrt(x));y=as.vector(sqrt(y))

  fit=loess(x~y, span = 0.99, family = "s")
  if(plot){
    plot(fit$x, fit$y)
    lines(sort(fit$x), fit$fitted[order(fit$x)])
  }
  return(as.vector(fit$residuals))
}

projectionSQ=function(allbranch){
  nv=projection(t(allbranch), method="AVE", returnNV = T)
  proj=resid(sqrt(allbranch), cbind(rep(1,length(nv)),nv))
  return(proj)
}



checkOrder=function(tree1, tree2, plot=F){
  both=intersect(tree1$tip.label, tree2$tip.label)

  tree1=unroot(pruneTree(tree1, both))
  tree2=unroot(pruneTree(tree2, both))

  tmpe1=as.data.frame(tree1$edge)
  tmpe1[match(1:length(both),tmpe1[,2]),2]=tree1$tip
  tmpe2=as.data.frame(tree2$edge)
  tmpe2[match(1:length(both),tmpe2[,2]),2]=tree2$tip

  map=matchNodes(tree1,tree2, method = "descendant")
  n=length(both)
  im=match(tree1$tip, tree2$tip)
  map=rbind(map, cbind(1:n, im))
  map=map[order(map[,1]),]

  edge1remap=tree1$edge


  edge1remap[,1]=map[edge1remap[,1],2]
  edge1remap[,2]=map[edge1remap[,2],2]
  if(plot){
    tmpd1=tree1$edge.length
    tmpd2=tree2$edge.length
    nn=character(length(tree1$edge))
    iim1=match(1:length(tree1$tip.label), tree1$edge[,2])
    nn[iim1]=tree1$tip.label
    par(mfrow=c(1,3))
    plot(tree1, use.edge.length = T); plot(tree2, use.edge.length = T)
    plotWtext(tmpd1, tmpd2, nn[])
  }
  return(all(edge1remap[,1]==tree2$edge[,1]) &&all(edge1remap[,2]==tree2$edge[,2]))

  #   #show(tmpd1)

  #  return(cbind(tmpd1, tmpd2))
}


getAllCladeEdges=function(tree, AncEdge){
  node=tree$edge[AncEdge,2]
  #get descendants
  iid=getDescendants(tree, node)
  #find their edges
  iim=match(iid, tree$edge[,2])
  iim
}



#if(T){
getNV=function(name1, name2, treesObj, residfun=residLN, plot=T){
  report=treesObj[["report"]]
  both=names(which(colSums(report[c(name1,name2),])==2))
  show(length(both))
  mastertree=pruneTree(treesObj[["master"]], both)
  allbranch=matrix(nrow=0, ncol=length(mastertree$edge.length))
  for ( i in 1:(length(treesObj)-3)){
    if(sum(is.na(match(both, treesObj[[i]]$tip.label)))==0){
      tmptree=(pruneTree(treesObj[[i]], both, mastertree))

      # show(c(length(tmptree$edge.length), ncol(allbranch)))
      allbranch=rbind(allbranch, tmptree$edge.length)
    }
  }
  nv=projection(t(allbranch), method="AVE", returnNV = T)
  mastertree$edge.length=nv
  # par(mfrow=c(1,3))
  res=correlateTrees(treesObj[[name1]], treesObj[[name2]], mastertree, residfun=residfun, plot=plot)
  res$nv=nv
  res$master=mastertree
  return(res)
}


getProjection=function(treesObj, tree1, tree2, maxT=treesObj$numTrees){
  both=intersect(tree1$tip.label, tree2$tip.label)
  tree1=unroot(pruneTree(tree1, both))
  tree2=unroot(pruneTree(tree2, both))
  allreport=treesObj$report[1:maxT,both]

  ss=rowSums(allreport)
  iiboth=which(ss==length(both))
  torm=setdiff(treesObj$masterTree$tip.label, both)
  allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
  for ( k in 1:length(iiboth)){
    tmptree=rescaleTree(unroot(drop.tip(treesObj$trees[[iiboth[k]]], torm)))
    allbranch[k, ]=tmptree$edge.length
  }
  allbranch
}

getProjectionPaths=function(treesObj, tree1, tree2, maxT=treesObj$numTrees){
  both=intersect(tree1$tip.label, tree2$tip.label)
  tree1=unroot(pruneTree(tree1, both))
  tree2=unroot(pruneTree(tree2, both))
  allreport=treesObj$report[1:maxT,both]
  ss=rowSums(allreport)
  iiboth=which(ss==length(both))
  allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
  ee=edgeIndexRelativeMaster(tree1, treesObj$masterTree)
  ii= match(namePaths(ee,T), colnames(treesObj$paths))
  allbranch=treesObj$paths[iiboth,ii]
  allbranch=scaleMat(allbranch)
  allbranch
}


correlateTreesAll=function(treesObj,  usePaths=F, useIndex=F,maxn=NULL, maxDo){
  maxT=treesObj$numTrees
  if (is.null(maxDo)){

    maxDo=maxT*(maxT-1)
  }
  corout=matrix(nrow=maxT, ncol=maxT)
  message("10")
  if(is.null(maxn)){
    maxn=treesObj$report%*%t(treesObj$report)
  }
  message("20")
  done=0
  # todo=length(maxn[upper.tri(maxn)]>=10)
  todo=100
  message("30")
  #corout[maxn<10]=0
  message(40)
  diag(corout)=1
  #corout[lower.tri(corout)]=0
  message("Starting loop")
  for (i in 1:(maxT-1)){
    for(j in (i+1):maxT){
      #  show(c(i,j))
      if (is.na(corout[i,j]) || maxn[i,j]<11){
        t0=as.double(Sys.time())
        tree1=treesObj$trees[[i]]
        tree2=treesObj$trees[[j]]

        bothIndex=which(colSums(treesObj$report[c(i, j),])==2)
        both=intersect(tree1$tip.label, tree2$tip.label)
        if(!useIndex){
          tree1=unroot(pruneTree(tree1, both))
          tree2=unroot(pruneTree(tree2, both))
        }
        allreport=treesObj$report[,bothIndex]
        ss=rowSums(allreport)
        iiboth=which(ss==length(bothIndex))
        #  torm=setdiff(treesObj$masterTree$tip.label, both)
        # allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))

        t1=as.double(Sys.time())
        message(paste("10 took", t1-t0))
        t0=t1
        if(! usePaths){
          torm=setdiff(treesObj$masterTree$tip.label, both)
          allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
          for ( k in 1:length(iiboth)){
            tmptree=rescaleTree(unroot(drop.tip(treesObj$trees[[iiboth[k]]], torm)))
            allbranch[k, ]=tmptree$edge.length
          }
        }
        else{
          if(!useIndex){
            message("Here")
            ee=edgeIndexRelativeMaster(tree1, treesObj$masterTree)
            ii= match(namePaths(ee,T), colnames(treesObj$paths))
            allbranch=treesObj$paths[iiboth,ii]
          }
          else{
            allbranch=getBranch(treesObj, bothIndex)
          }
          t1=as.double(Sys.time())
          message(paste("20 took", t1-t0))
          t0=t1
          allbranch=scaleMat(allbranch)
        }

        message("done")
        nb=length(both)
        proj=t(projection(t(allbranch), method="AVE", returnNV = F))
        # i1=match(i, iiboth)
        #j1=match(j,iiboth)
        #corout[i,j]=cor(proj[i1, ], proj[j1,])
        #  tmpcor=cor(t(proj))

        ai=which(maxn[iiboth, iiboth]==nb, arr.ind = T)
        t1=as.double(Sys.time())
        message(paste("30 took", t1-t0))
        t0=t1
        for (m in 1:nrow(ai)){
          k=sort(ai[m,])[1]
          l=sort(ai[m,])[2]

          tmpcor=cor(proj[k,], proj[l,])
          if (is.na(tmpcor)){
            tmpcor=0
          }
          corout[iiboth[k], iiboth[l]]=tmpcor


        }
        t1=as.double(Sys.time())
        message(paste("40 took", t1-t0))
        t0=t1
        done=done+nrow(ai)
        message(paste("Done with",done, "out of", todo))
        #  message(paste(sum(is.na(corout))," left"), appendLF = T)
        if(done>=maxDo){
          message("DOne")
          return(corout)
        }
        #generate the projection
      }
    }
  }

}



#assume the binTree is already in canonical form
correlateTreesBinary=function(treesObj,  binTree, usePaths=F, maxDo=NULL, species.list=NULL, useSQ=F){
  maxT=treesObj$numTrees
  if (is.null(maxDo)){

    maxDo=maxT
  }
  corout=matrix(nrow=maxT, ncol=1)
  pout=matrix(nrow=maxT, ncol=1)
  rownames(corout)=rownames(pout)=names(treesObj$trees)
  show(binTree$tip.label )
  binReport=as.vector(as.numeric(binTree$tip.label %in% colnames(treesObj$report)))
  show((binReport))
  names(binReport)=colnames(treesObj$report)
  maxn=treesObj$report[, species.list]%*%(binReport[species.list])

  done=0
  todo=length(maxn>=10)
  corout[maxn<10]=0

  for (i in 1:maxT){

    if (is.na(corout[i,1])){
      tree1=treesObj$trees[[i]]
      if(! is.null(species.list)){
        tree1=unroot(pruneTree(tree1, species.list))
      }
      both=tree1$tip.label
      bothIndex=match(both, colnames(treesObj$report))
      allreport=treesObj$report[,bothIndex]
      ss=rowSums(allreport)
      iiboth=which(ss==length(both))
      #  torm=setdiff(treesObj$masterTree$tip.label, both)
      binTreeUse=unroot(pruneTree(binTree,tree1$tip.label))
      allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
      if(length(both)<10){
        next
      }

      if(! usePaths){
        torm=setdiff(treesObj$masterTree$tip.label, both)
        allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
        for ( k in 1:length(iiboth)){
          tmptree=rescaleTree(unroot(drop.tip(treesObj$trees[[iiboth[k]]], torm)))
          allbranch[k, ]=tmptree$edge.length
        }
      }
      else{

        ii= match(namePaths(edgeIndexRelativeMaster(tree1, treesObj$masterTree),T), colnames(treesObj$paths))
        ii2=match(namePaths(edgeIndexRelativeMaster(binTreeUse, treesObj$masterTree),T), colnames(treesObj$paths))
        #show(ii)
        #show(ii2)
        stopifnot(all(ii=ii2))
        allbranch=treesObj$paths[iiboth,ii]

        allbranch=scaleMat(allbranch)
      }
      # message("done")
      nb=length(both)
      if(!useSQ){
        proj=t(projection(t(allbranch), method="AVE", returnNV = F))
      }
      else{
        proj=projectionSQ(allbranch)
      }
      # i1=match(i, iiboth)
      #j1=match(j,iiboth)
      #corout[i,j]=cor(proj[i1, ], proj[j1,])
      #  tmpcor=cor(t(proj))
      ai=which(maxn[iiboth, 1]==nb)
      # show(iiboth[1])


      tmp=simpleAUCmat(binTreeUse$edge.length, (proj[ai, ,drop=F]))

      corout[iiboth[ai]]=tmp$auc

      pout[iiboth[ai]]=tmp$pp
      done=done+length(ai)
      show(length(ai))
      #  message(paste("Done with",done, "out of", todo))
      #  message(paste(sum(is.na(corout))," left"), appendLF = T)
      if(done>=maxDo){
        message("DONE")
        return(list(r=corout, p=pout))
      }
      #generate the projection
    }
  }

  return(list(r=corout, p=pout))
}


plotTreesBinary=function(treesObj,  binTree, index, species.list=NULL){
  maxT=treesObj$numTrees

  binReport=as.vector(as.numeric(binTree$tip.label %in% colnames(treesObj$report)))

  names(binReport)=colnames(treesObj$report)
  maxn=treesObj$report[, species.list]%*%(binReport[species.list])




  tree1=treesObj$trees[[index]]
  if(! is.null(species.list)){
    tree1=unroot(pruneTree(tree1, species.list))
  }
  both=tree1$tip.label
  bothIndex=match(both, colnames(treesObj$report))
  allreport=treesObj$report[,bothIndex]
  ss=rowSums(allreport)
  iiboth=which(ss==length(both))
  #  torm=setdiff(treesObj$masterTree$tip.label, both)
  binTreeUse=unroot(pruneTree(binTree,tree1$tip.label))
  allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))

  ii= match(namePaths(edgeIndexRelativeMaster(tree1, treesObj$masterTree),T), colnames(treesObj$paths))
  ii2=match(namePaths(edgeIndexRelativeMaster(binTreeUse, treesObj$masterTree),T), colnames(treesObj$paths))
  plot(tree1,use.edge.length = F)
  plot(binTreeUse,use.edge.length = F)
  show(cbind(ii,ii2))
  stopifnot(all(ii==ii2))
  allbranch=treesObj$paths[iiboth,ii]
  thisgene=which(iiboth==index)
  allbranch=scaleMat(allbranch)

  nb=length(both)

  proj=t(projection(t(allbranch), method="AVE", returnNV = F))
  nv=t(projection(t(allbranch), method="AVE", returnNV = T))

  plot(nv,proj[thisgene,], col=binTreeUse$edge.length+1)

}



correlateTrees=function(tree1, tree2, mastertree, residfun=residLN, plot=F, cutoff=0.00001, Tree1Bin=F){
  both=intersect(tree1$tip.label, tree2$tip.label)
  if(length(both)<10){
    return(0)
  }
  iibad1=which(tree1$edge.length<cutoff)
  iibad2=which(tree2$edge.length<cutoff)
  show(c(length(iibad1), length(iibad2)))
  show(tree1$edge.length)
  if (!Tree1Bin){
    tree1=rescaleTree(tree1)
  }
  tree2=rescaleTree(tree2)
  tree1$edge.length[iibad1]=mastertree$edge.length[iibad1]
  tree2$edge.length[iibad2]=mastertree$edge.length[iibad2]
  if(!Tree1Bin){
    e1=residfun(t(tree1$edge.length), mastertree$edge.length, plot=F)
  }
  else{
    e1=tree1$edge.length
  }
  e2=residfun(t(tree2$edge.length), mastertree$edge.length, plot=F)
  cc=cor((e1), (e2))
  nn=character(length(e1))
  iim=match(1:length(tree1$tip.label), tree1$edge[,2])

  nn[iim]=tree1$tip.label

  if(plot){


    plotWtext(e1, e2, nn)
    title(paste("R=", round(cc,2)))
  }

  return(list(l1=tree1$edge.length, l2=tree2$edge.length,e1=e1, e2=e2, cor=cc, names=nn, tree1=tree1, tree2=tree2))

}



correlateTreesProj=function(treeIn1, treeIn2, treesObj, residfun=residLN, plot=F, cutoff=-1, usePaths=T, tree1Bin=F, useIndex=F, species.list=NULL){
  if(is.character(treeIn1)){
    tree1=treesObj$trees[[treeIn1]]
  }
  else{
    tree1=treeIn1
  }
  if(is.character(treeIn2)){
    tree2=treesObj$trees[[treeIn2]]
  }
  else{
    tree2=treeIn2
  }
  both=intersect(tree1$tip.label, tree2$tip.label)
  if(!is.null(species.list)){
    both=intersect(both, species.list)
  }


  torm=setdiff(treesObj$masterTree$tip.label, both)
  tree1=pruneTree(tree1, both)
  tree1=unroot(tree1)
  if(tree1Bin){ #fix any edges that were created through pruning
    tree1$edge.length[tree1$edge.length>1]=1
  }
  tree2=pruneTree(tree2, both)
  tree2=unroot(tree2)
  allreport=treesObj$report[,both]
  ss=rowSums(allreport)
  iiboth=which(ss==length(both))
  if (! usePaths){
    allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
    for ( i in 1:length(iiboth)){
      tmptree=rescaleTree(drop.tip(treesObj$trees[[iiboth[i]]], torm))
      allbranch[i, ]=tmptree$edge.length
    }

  }
  else{
    if(! useIndex){
      ee=edgeIndexRelativeMaster(tree1, treesObj$masterTree)
      ii= match(namePaths(ee,T), colnames(treesObj$paths))
      allbranch=treesObj$paths[iiboth,ii]
      show(sum(is.na(allbranch)))
      allbranch=scaleMat(allbranch)

      nv=projection(t(allbranch), method="AVE", returnNV = T)
      mastertree=treesObj$master
      mastertree$edge.length=nv
      res=correlateTrees(tree1, tree2, mastertree, residfun=residfun, plot=plot, cutoff=cutoff, Tree1Bin=tree1Bin)

      res$nv=nv
      res$allbranch=allbranch
    }
    else{
      allbranch=getBranch(treesObj, bothIndex)
      show(rownames(allbranch)[1])
      allbranch=scaleMat(allbranch)
      nv=projection(t(allbranch), method="AVE",returnNV = T)
      rr=resid(allbranch, model.matrix(~0+nv))
      rownames(rr)=rownames(allbranch)
      show(dim(rr))
      show(rownames(allbranch)[1])
      plot(rr[name1,], rr[name2,])
      res=list()
    }
  }

  return(res)
}



correlateTreesAll1=function(name1, name2, treesObj, residfun=residLN, plot=F, cutoff=-1, usePaths=F){
  tree1=treesObj$trees[[name1]]
  tree2=treesObj$trees[[name2]]
  both=intersect(tree1$tip.label, tree2$tip.label)
  torm=setdiff(treesObj$mastertree$tip.lables, both)
  tree1=pruneTree(tree1, both)
  tree2=pruneTree(tree2, both)
  allreport=treesObj$report[,both]
  ss=rowSums(allreport)

  iiboth=which(ss==length(both))
  if (! usePaths){
    allbranch=matrix(nrow=length(iiboth), ncol=length(tree1$edge.length))
    for ( i in 1:length(iiboth)){
      # tmptree=rescaleTree(pruneTree(treesObj$trees[[iiboth[i]]], both))
      tmptree=rescaleTree(drop.tip(treesObj$trees[[iiboth[i]]], torm))
      #  show(c(length(tmptree$edge.lenght), ncol(allbranch)))
      allbranch[i, ]=tmptree$edge.length
    }
  }
  else{
    allbranch=treesObj$paths[iiboth,getEdgeIndex(tree1, treesObj$masterTree)]
    #allbranch=scaleMat(allbranch)
  }
  #nv=projection(t(allbranch), method="AVE", returnNV = T)
  nv=colMeans(allbranch)
  nn=names(iiboth)
  rownames(allbranch)=names(iiboth)
  allbranch=resid(allbranch, model.matrix(~1+nv))
  cc=cor(t(allbranch))
  show(dim(cc))
  show(length(nn))
  ii=which(trees$inter[nn,nn]==length(both))
  return(list(cc,ii))
}




mapEdge=function(tree1, tree2){
  map=matchNodes(tree1,tree2, method = "descendant")
  n=length(tree1$tip)
  im=match(tree1$tip, tree2$tip)
  map=rbind(map, cbind(1:n, im))
  map=map[order(map[,1]),]
  #show(n)
  #show(map)
  edge1remap=tree1$edge
  #show(nrow(edge1remap))
  #show(nrow(tree2$edge))

  edge1remap[,1]=map[edge1remap[,1],2]
  edge1remap[,2]=map[edge1remap[,2],2]
  edge1remap

}


plotContinuousCharXY=function(gene, treesObj, tip.vals, tip.vals.ref=NULL,  col=NULL, residfun=residLO, useDiff=T, xlab){
  #get the tree projection
  tip.vals=tip.vals[!is.na(tip.vals)]

  stopifnot(gene %in% names(treesObj$trees))
  tree=treesObj$trees[[gene]]
  stopifnot(!is.null(names(tip.vals)))
  both=intersect(tree$tip.label, names(tip.vals))

  stopifnot(length(both)>10)


  torm=setdiff(treesObj$masterTree$tip.label, both)
  tree=pruneTree(tree, both)
  tip.vals=tip.vals[both]
  allreport=treesObj$report[,both]
  ss=rowSums(allreport)
  iiboth=which(ss==length(both))


  ee=edgeIndexRelativeMaster(tree, treesObj$masterTree)
  ii= match(namePaths(ee,T), colnames(treesObj$paths))

  allbranch=treesObj$paths[iiboth,ii]

  allbranch=scaleMat(allbranch)
  show(sum(is.na(allbranch)))
  nv=projection(t(allbranch), method="AVE", returnNV = T)

  proj=residfun(tree$edge.length, nv)
  show(length(tree$edge.length))
  treeChar=edgeVarsDiff(tree, tip.vals)
  show(length(treeChar))
  show(length(proj))
  nn=nameEdges(tree)
  nn[nn!=""]=speciesNames[nn[nn!=""], ]
  #par(mfrow=c(2,2), mai=rep(0.7,4))
  #plotWtext(sqrt(nv), sqrt(tree$edge.length), xlab="char", ylab="Gene branch length", labels = nn)

  plotWtext(treeChar$edge.length, proj, xlab=xlab, ylab="relative gene branch length", labels = nn)
  stat=cor.test(treeChar$edge.length, proj, method="s")
  mtext(gene,side = 3, line=2, font=2)
  mtext(paste0("r=", round(stat$estimate,2), ";  p-value=", format.pval(stat$p.value)), side = 3, line=0.5, cex=.7)

  if(!is.null(tip.vals.ref)){
    treeCharRef=edgeVars(tree, tip.vals.ref, useDiff=useDiff)
    proj=resid(rbind(proj), model.matrix(~1+treeCharRef$edge.length))[1,]
  }




}

plotWtext=function(x,y, labels, text.cex=0.7, ...){plot(x,y, pch=19, col="#00008844", xlim=range(x)+c(0,0.7),...); textplot(x,y, words=labels, cex=text.cex)}


plotResidsVsChar=function(x,y, labels, text.cex=0.7, names=T,...){
  ii=which(!is.na(x)&!is.na(y))
  show(range(x[ii]))
  nn=names(x)[ii]
  nn=speciesNames[nn,1]
  nn[is.na(nn)]=""

  plot(x[ii], y[ii], col="#0000AAAA",xlab="RER", ...);
  if(names)
    textplot(x[ii],y[ii], words=nn, cex=text.cex,new = F)
}


#' Calculates Rho-signed negative log-base-ten p-value for use in enrichment functions

#' @param res The output from RERconverge correlation functions (correlateWithContinuousPhenotype, correlateWithBinaryPhenotype, getAllCor)
#' @return A dataframe of Rho-signed negative log-base-ten p-values for all genes, NAs removed
#' @export

getStat=function(res){
  stat=sign(res$Rho)*(-log10(res$P))
  names(stat)=rownames(res)
  #deal with duplicated genes
  genenames=sub("\\..*", "",names(stat))
  multname=names(which(table(genenames)>1))
  for(n in multname){
    ii=which(genenames==n)
    iimax=which(max(stat[ii])==max(abs(stat[ii])))
    stat[ii[-iimax]]=NA
  }
  sum(is.na(stat))
  stat=stat[!is.na(stat)]

  stat
}

varExplainedWithNA=function (dat, val, adjust=T)
{

  adje=vare=double(nrow(dat))
  names(vare)=rownames(dat)
  for(i in 1:nrow(dat)){
    ii=which(!is.na(dat[i,])&!is.na(val))
    if(length(ii)>10){
      mod0 = cbind(rep(1, length(ii)))
      mod=model.matrix(~1+val[ii])
      n=length(ii)

      adj=(n-1)/(n-ncol(mod))

      evince
      resid = resid(dat[i,ii, drop=F], mod)
      resid0 = resid(dat[i,ii,drop=F], mod0)
      rss1 = resid^2 %*% rep(1, n)
      rss0 = resid0^2 %*% rep(1, n)

      vare[i]=1-rss1/rss0*adj
      adje[i]=adj        }
  }
  return(cbind(vare, adje))
}








findJoining=function(matIndex, speciesIndex){
  ss=colSums(matIndex[speciesIndex,])
  ssi=which(ss==2)
  ii=ssi[which(colSums(matIndex[ssi,ssi, drop=F])==0)]
  res=matrix(nrow=length(ii), ncol=5)
  count=0
  for (i in ii){
    tmpi=which(matIndex[speciesIndex,i]==1)
    count=count+1
    res[count,]=c(i, speciesIndex[tmpi],tmpi)
  }
  res
}
findAllPaths=function(matIndex, speciesIndex, key.species=1, tri.node=NULL){
  if(is.null(tri.node)){
    tri.node=which(colSums(matIndex)%%2==1)
  }
  speciesIndex=speciesIndex[-which(speciesIndex==key.species)]
  count=0
  paths=matrix(ncol=2, nrow=0)
  while(length(speciesIndex)>1 &&count<200){
    count=count+1
    ss=colSums(matIndex[speciesIndex,])
    res=findJoining(matIndex, speciesIndex)
    toRm=integer()
    for ( i in 1:nrow(res)){
      ti=speciesIndex[tmpti<-which(matIndex[speciesIndex, i]==1)]

      paths=rbind(paths, c(res[i,1], res[i,2]))
      paths=rbind(paths, c(res[i,1], res[i,3]))
    }
    toRm=as.vector(res[, 4:5])
    speciesIndex=speciesIndex[-toRm]
    speciesIndex=c(speciesIndex, res[,1])
  }
  paths=rbind(paths, c(tri.node, key.species))
  paths
}

getBranch=function(treesObj, speciesIndex,key.species=1, tri.node=NULL){
  res=findAllPaths(treesObj$matAnc, speciesIndex)
  ii=which(rowSums(treesObj$report[,speciesIndex])>=length(speciesIndex))
  allBranch=treesObj$paths[ii,treesObj$matIndex[res[,c(2:1)]]]

}









getAncestor=function(tree, nodeN){
  if(is.character(nodeN)){
    nodeN=which(tree$tip.label==nodeN)
  }
  im=which(tree$edge[,2]==nodeN)
  return(tree$edge[im,1])
}






matchNodesInjectUpdate=function (tr1, tr2){
  ancMatFrom=getAncestorMatrix(tr1)


  ancMatTo=getAncestorMatrix(tr2) # this is the

  outMat=matrix(0,nrow=nrow(ancMatFrom), ncol=ncol(ancMatTo))
  colnames(outMat)=colnames(ancMatTo)
  #rownames(outMat)=rownames(ancMatFrom)
  outMat[, colnames(ancMatFrom)]=ancMatFrom
  tt=tcrossprod(outMat,ancMatTo)
  tt2=sqrt(outer(rowSums(outMat), rowSums(ancMatTo[, colnames(ancMatFrom)])))
  ii=which(tt/tt2==1, arr.ind = T)

  rr2<-rowSums(ancMatTo)
  #ii is mappint tr1 to tr2
  iidouble=which(table(ii[,1])>1)
  if(length(iidouble)>0){
    for (i in iidouble){
      iiresolve=which(ii[,1]==i)

      #     show(rr2[ii[iiresolve,2]])
      #  show(ii[iiresolve,])
      j=which.min(rr2[ii[iiresolve,2]])
      #show(i)
      j=ii[iiresolve[j],2]
      #show(j)
      iirm=which(ii[,1]==i&ii[,2]!=j)
      #show(ii[iirm,])

      ii[iirm,]=NA
      #show(ii[iiresolve,])
    }
  }
  ii[,1]=ii[,1]+length(tr1$tip.label)
  ii[,2]=ii[,2]+length(tr2$tip.label)

  ii=ii[!is.na(ii[,1]),]

  ii=rbind(cbind(1:length(tr1$tip.label), match(tr1$tip.label, tr2$tip.label)),ii)
  if(nrow(ii)!=length(tr1$tip.label)+tr1$Nnode){
    stop("Discordant tree topology detected")
  }
  ii
  #ii=ii[order(ii[,1]),]
}

