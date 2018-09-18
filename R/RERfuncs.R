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

#' reads trees from a 2 column , tab seperated, file
#' The first columns is the gene name and the second column is the corresponding tree in parenthetic format known as the Newick or New Hampshire format

#' @param file The path to the tree file
#' @param  max.read this function takes a while for  a whole genome so max.read is useful for testing
#' @return A trees object of class "treeObj"
#' @export
readTrees=function(file, max.read=NA){
  tmp=scan(file, sep="\t", what="character")
  trees=vector(mode = "list", length = min(length(tmp)/2,max.read, na.rm = T))
  treenames=character()
  maxsp=0; # maximum number of species

  for ( i in 1:min(length(tmp),max.read*2, na.rm = T)){
    if (i %% 2==1){
      treenames=c(treenames, tmp[i])
    }
    else{
      trees[[i/2]]=unroot(read.tree(text=tmp[i]))
      #check if it has more species
      if(length(trees[[i/2]]$tip.label)>maxsp){
        maxsp=length(trees[[i/2]]$tip.label)
        allnames=trees[[i/2]]$tip.label
      }
    }

  }
  names(trees)=treenames
  treesObj=vector(mode = "list")
  treesObj$trees=trees
  treesObj$numTrees=length(trees)
  treesObj$maxSp=maxsp

  message(paste("max is ", maxsp))

  report=matrix(nrow=treesObj$numTrees, ncol=maxsp)
  colnames(report)=allnames

  rownames(report)=treenames
  for ( i in 1:nrow(report)){
    ii=match(allnames, trees[[i]]$tip.label)
    report[i,]=1-is.na(ii)

  }
  treesObj$report=report



  ii=which(rowSums(report)==maxsp)

  #Create a master tree with no edge lengths
  master=trees[[ii[1]]]
  master$edge.length[]=1
  treesObj$masterTree=master




  treesObj$masterTree=rotateConstr(treesObj$masterTree, sort(treesObj$masterTree$tip.label))
  #this gets the abolute alphabetically constrained order when all branches
  #are present
  tiporder=treeTraverse(treesObj$masterTree)

  #treesObj$masterTree=CanonicalForm(treesObj$masterTree)

  for ( i in 1:treesObj$numTrees){

    treesObj$trees[[i]]=rotateConstr(treesObj$trees[[i]], tiporder)

  }



  ap=allPaths(master)
  treesObj$ap=ap
  matAnc=(ap$matIndex>0)+1-1
  matAnc[is.na(matAnc)]=0

  paths=matrix(nrow=treesObj$numTrees, ncol=length(ap$dist))
  for( i in 1:treesObj$numTrees){
    paths[i,]=allPathMasterRelative(treesObj$trees[[i]], master, ap)
  }
  paths=paths+min(paths[paths>0], na.rm=T)
  treesObj$paths=paths
  treesObj$matAnc=matAnc
  treesObj$matIndex=ap$matIndex
  treesObj$lengths=unlist(lapply(treesObj$trees, function(x){sqrt(sum(x$edge.length^2))}))

  ii=which(rowSums(report)==maxsp)
  if(length(ii)>20){
    message (paste0("estimating master tree branch lengths from ", length(ii), " genes"))
    tmp=lapply( treesObj$trees[ii], function(x){x$edge.length})

    allEdge=matrix(unlist(tmp), ncol=2*maxsp-3, byrow = T)
    allEdge=scaleMat(allEdge)
    allEdgeM=apply(allEdge,2,mean)
    treesObj$masterTree$edge.length=allEdgeM
  }
  else{
    message("Not enough genes with all species present: master tree has no edge.lengths")
  }
  colnames(treesObj$paths)=namePathsWSpecies(treesObj$masterTree)
  class(treesObj)=append(class(treesObj), "treesObj")
  treesObj
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
    breaks=round(breaks,3)
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
    stop("Discordant tree topology detected")
  }

  Nodes
}



#' @keywords  internal
allPaths=function(tree){
  dd=dist.nodes(tree)
  allD=double()
  nn=matrix(nrow=0, ncol=2)
  nA=length(tree$tip.label)+tree$Nnode
  matIndex=matrix(nrow=nA, ncol=nA)
  index=1
  for ( i in 1:nA){
    ia=getAncestors(tree,i)
    if(length(ia)>0){
      allD=c(allD, dd[i, ia])
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
#' @param method Method used to compute correlations. Accepts the same arguments as \code{\link{cor}}. Set to "auto" to select automatically based on the number of unique values in charP. This will also auto set the winsorization for Pearson correlation. Set winsorize=some number to override
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
#' @param method Method used to compute correlations. Accepts the same arguments as \code{\link{cor}}. Set to "auto" to select automatically based on the number of unique values in charP. This will also auto set the winsorization for Pearson correlation. Set winsorize=some number to override
#' @param min.sp Minimum number of species that must be present for a gene
#' @param winsorize Winsorize values before computing Pearson correlation. Winsorize=3, will set the 3 most extreme values at each end to the the value closest to 0.
#' @export

correlateWithContinuousPhenotype=function(RERmat,charP, min.sp=10,  winsorize=3){
  getAllCor(RERmat, charP, min.sp, min.pos=0, method = "p", winsorize = winsorize )
}




#'Computes the association statistics between RER from \code{\link{getAllResiduals}} and a phenotype paths vector made with \code{\link{tree2Paths}} or \code{\link{char2Paths}}
#' @param RERmat RER matrix returned by \code{\link{getAllResiduals}}
#' @param charP phenotype vector returned by \code{\link{tree2Paths}} or \code{\link{char2Paths}}
#' @param method Method used to compute correlations. Accepts the same arguments as \code{\link{cor}}. Set to "auto" to select automatically based on the number of unique values in charP. This will also auto set the winsorization for Pearson correlation. Set winsorize=some number to override
#' @param min.sp Minimum number of species that must be present for a gene
#' @param min.pos Minimum number of species that must be present in the foreground (non-zero phenotype values)
#' @param winsorize Winsorize values before computing Pearson correlation. Winsorize=3, will set the 3 most extreme values at each end to the the value closest to 0.
#' @param weighted perform weighted correlation. This option needs to be set if the clade weights computed in \code{\link{foreground2Tree}(wholeClade=T)} are to be used. This setting will treat the clade a single observation for the purpose of p-value estimation.
#' @note  winsorize is in terms of number of observations at each end, NOT quantiles
#' @return A list object with correlation values, p-values, and the number of data points used for each tree
#' @export
getAllCor=function(RERmat, charP, method="auto",min.sp=10, min.pos=2, winsorize=NULL,weighted=F){
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
      if(is.null(winsorize)){
        message("Setting winsorise=3")
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

  colnames(corout)=c("Rho", "N", "P")
  if(!is.null(winsorize)){
    charP=win(charP, winsorize)
  }

  for( i in 1:nrow(corout)){

    if(((nb<-sum(ii<-(!is.na(charP)&!is.na(RERmat[i,]))))>=min.sp)){
      if (method!="p"&&sum(charP[ii]!=0)<min.pos){
        next
      }

      if(!weighted){


        if (!is.null(winsorize)){
          x=win(RERmat[i,], winsorize)
        }
        else{
          x=RERmat[i,]
        }
        cres=cor.test(x, charP, method=method, exact=F)
        corout[i,1:3]=c(cres$estimate, nb, cres$p.value)
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
  corout
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
#' @return A numer of trees by number of paths matrix of relative evolutionary rates. Only an independent set of paths has non-NA values for each tree.
#' @export
getAllResiduals=function(treesObj, cutoff=NULL, transform="sqrt", weighted=T,  useSpecies=NULL,  min.sp=10, scale=T,  doOnly=NULL, maxT=NULL, scaleForPproj=F, mean.trim=0.05){

  if(is.null(cutoff)){
    cutoff=quantile(treesObj$paths, 0.05, na.rm=T)
    message(paste("cutoff is set to", cutoff))
  }
  if (weighted){
    weights=computeWeightsAllVar(treesObj$paths, transform=transform, plot=T)
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
      iiboth=which(ss==length(both))

      nb=length(both)
      ai=which(maxn[iiboth]==nb)


      message(paste("i=", i))


      if(T){


        ee=edgeIndexRelativeMaster(tree1, treesObj$masterTree)

        ii= treesObj$matIndex[ee[, c(2,1)]]

        allbranch=treesObj$paths[iiboth,ii]
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
foreground2Paths = function(foreground,treesObj, plotTree=F, clade=c("ancestral","terminal","all","weighted"), useSpecies=NULL){
  #res = treesObj$masterTree
  #res$edge.length <- rep(0,length(res$edge.length))
  #res$edge.length[nameEdges(treesObj$masterTree) %in% foreground] = 1
  #names(res$edge.length) = nameEdges(treesObj$masterTree)
  res = foreground2Tree(foreground, treesObj, plotTree=plotTree, clade=clade, useSpecies=useSpecies)
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
#' of the strings "ancestral", "terminal", "all", or "weighted".
#' @param useSpecies Give only a subset of the species to use for ancestral state reconstruction
#' (e.g., only those species for which the trait can be reliably determined).
#' @return A tree with edge.lengths representing phenotypic states
#' @export
foreground2Tree = function(foreground,treesObj, collapse2anc=T, plotTree=T,  wholeClade=F, clade=c("ancestral","terminal","all","weighted"), useSpecies=NULL){
  clade <- match.arg(clade) #should error if not an allowed option
  wholeClade = T
  collapse2anc = T
  if (clade %in% c("ancestral","terminal")) {
    wholeClade = F
  }
  #if(wholeClade){
  #  collapse2anc=T
  #}
  if (clade == "terminal") {
    collapse2anc = F
  }
  res = treesObj$masterTree
  if (!is.null(useSpecies)) {
    useSpecies = intersect(useSpecies, res$tip.label)
    res = pruneTree(res, useSpecies)
    sp.miss = setdiff(res$tip.label, useSpecies)
    if (length(sp.miss) > 0) {
      message(paste0("Species from useSpecies not present in master tree: ", paste(sp.miss,
                                                                                   collapse = ",")))
    }
  } else {
    useSpecies = res$tip.label
  }
  foreground = intersect(foreground, useSpecies)
  res$edge.length <- rep(0,length(res$edge.length))
  if(!collapse2anc){
    res$edge.length[nameEdges(treesObj$masterTree) %in% foreground] = 1
    names(res$edge.length) = nameEdges(treesObj$masterTree)
  }
  else{
    tip.vals=rep(0, length(treesObj$masterTree$tip.label))
    names(tip.vals)=treesObj$masterTree$tip.label
    tip.vals[foreground]=1
    tmp=cbind(as.character(tip.vals))
    rownames(tmp)=names(tip.vals)
    tip.vals=tmp

    ancres=ancestral.pars(res, df<-as.phyDat(tip.vals, type="USER", levels=unique(as.character(tip.vals))),type="ACCTRAN" )

    ancres=unlist(lapply(ancres, function(x){x[2]}))
    internalVals=ancres
    evals=matrix(nrow=nrow(treesObj$masterTree$edge), ncol=2)
    eres=ancres
    evals[,1]=eres[treesObj$masterTree$edge[,1]]
    evals[,2]=eres[treesObj$masterTree$edge[,2]]
    res$edge.length=evals[,2]-evals[,1]

    res$edge.length[res$edge.length<1]=0
  }
  if(wholeClade){ #should implement for "all" and "weighted"
    edgeIndex=which(res$edge.length>0)
    for(i in edgeIndex){
      clade.edges=getAllCladeEdges(res, i)
      clade.edges=unique(c(i, clade.edges))
      if (clade == "weighted") {
        res$edge.length[clade.edges]=1/length(clade.edges)
      } else {
        res$edge.length[clade.edges]=1
      }
    }
  }
  if(plotTree){
    res2=res
    mm=min(res2$edge.length[res2$edge.length>0])
    res2$edge.length[res2$edge.length==0]=max(0.02,mm/20)
    plot(res2)
  }
  res
}




#' @keywords internal
nameEdges=function(tree){
  nn=character(nrow(tree$edge))
  iim=match(1:length(tree$tip.label), tree$edge[,2])
  nn[iim]=tree$tip.label
  nn
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
tree2Paths=function(tree, treesObj, binarize=NULL, useSpecies=NULL){
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

  treePaths=allPaths(tree)
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



if(F){
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
