combinePermData=function(permdat1, permdat2){
  #combine results
  allres=merge(permdat1[[1]],permdat2[[1]],by="row.names")
  rownames(allres)=allres[,1]
  allres[,1]=NULL

  #combine enrich
  enrich1=permdat1[[2]]
  enrich2=permdat2[[2]]
  groups=length(enrich1)
  allnewenrich=vector("list", groups)
  count=1
  while(count<=groups){
    newenrich=merge(enrich1[[count]], enrich2[[count]], by="row.names")
    rownames(newenrich)=newenrich[,1]
    newenrich[,1]=NULL
    allnewenrich[[count]]=newenrich
    count=count+1
  }
  data=vector("list", 2)
  data[[1]]=allres
  data[[2]]=allnewenrich
  data
}




#'Calculates permuted correlation and enrichment statistics
#' @param numperms An integer number of permulations
#' @param traitvec A named phenotype vector
#' @param RERmat An RER matrix calculated using \code{\link{getAllResiduals}}
#' @param annotlist Pathway annotations
#' @param treetop A rooted, fully dichotomous tree derived from the treesObj master tree from \code{\link{readTrees}}.  Must not contain species not in traitvec
#' @param trees treesObj from \code{\link{readTrees}}
#' @param type One of "simperm", "sim", or "perm" for permulations, simulations, or permutations, respectively
#' @param winR Integer winzorization value for RERs
#' @param winT Integer winzorization value for trait
#' @note  winsorize is in terms of number of observations at each end, NOT quantiles
#' @return A list object with enrichment statistics, correlation p-val, rho, and correlation effect size
#' @export
getPerms=function(numperms, traitvec, RERmat, annotlist, treetop, trees, type="simperm", winR=NULL, winT=NULL){
  #get real enrich and cors
  vec=traitvec
  realpaths=RERconverge::char2Paths(vec, trees)
  realresults=RERconverge::getAllCor(RERmat, realpaths, method = "p", min.pos = 0, winsorizeRER = winR, winsorizetrait=winT)
  realstat=getStat(realresults)
  realenrich=fastwilcoxGMTall(realstat, annotlist, outputGeneVals=T)

  #sort real enrichments
  groups=length(realenrich)
  c=1
  while(c<=groups){
    current=realenrich[[c]]
    realenrich[[c]]=current[order(rownames(current)),]
    c=c+1
  }

  #make enrich list/matrices to fill
  permutationresults=data.frame(matrix(ncol=numperms, nrow=nrow(realresults)))
  rownames(permutationresults)=rownames(realresults)
  rhos=data.frame(matrix(ncol=numperms, nrow=nrow(realresults)))
  rownames(rhos)=rownames(realresults)
  permutationenrichments=vector("list", length(realenrich))
  c=1
  while(c<=length(realenrich)){
    newdf=data.frame(matrix(ncol=numperms, nrow=nrow(realenrich[[c]])))
    rownames(newdf)=rownames(realenrich[[c]])
    permutationenrichments[[c]]=newdf
    c=c+1
  }

  v=vec
  rr=RERmat
  meredithplustreesv2=trees

  #pick number of rows of stat matrix
  st=getStat(realresults)
  r=length(st)
  signpvals=data.frame(matrix(ncol=numperms, nrow=r))
  rownames(signpvals)=names(st)

  counter=1;
  while(counter<=numperms){

    print(paste("running permutation: ", counter))

    #get correlation results
    out=getNullCor(v, rr, treetop, meredithplustreesv2, type = type, winR=winR, winT=winT)
    permutationresults[,counter]=out$P
    rhos[,counter]=out$Rho

    #get enrich results
    stat=getStat(out)
    #store stat values
    signpvals[,counter]=stat

    enrich=fastwilcoxGMTall(stat, annotlist, outputGeneVals=T)

    #sort and store enrichment results
    groups=length(enrich)
    c=1
    while(c<=groups){
      current=enrich[[c]]
      enrich[[c]]=current[order(rownames(current)),]
      permutationenrichments[[c]][,counter]=enrich[[c]]$stat
      c=c+1
    }

    counter=counter+1
  }
  data=vector("list", 4)
  data[[1]]=permutationresults
  data[[2]]=permutationenrichments
  data[[3]]=rhos
  data[[4]]=signpvals
  names(data)=c("correlationpvals", "signedenrichmentstats", "rhos", "correlationeffectsize")
  data
}


#'Calculates permutation pvals from output of \code{\link{getPerms}}
#' @param realenrich Real enrichment statistics from \code{\link{fastwilcoxGMTall}}
#' @param permenrich signedenrichmentstats from \code{\link{getPerms}}
#' @return A list object with permulation p-values
#' @export
permpvalenrich=function(realenrich, permenrich){
  #takes real and perm enrich with rows in the same order
  enrichpvals=vector("list", length(realenrich))
  groups=length(realenrich)
  count=1
  while(count<=groups){
    currreal=realenrich[[count]]
    currenrich=permenrich[[count]]
    rowlen=nrow(currenrich)
    rowcount=1
    pvallist=c()
    while(rowcount<=rowlen){
      lessnum=sum(abs(currenrich[rowcount,])>abs(currreal[rowcount,]$stat))
      pval=lessnum/ncol(currenrich)
      pvallist=c(pvallist, pval)
      rowcount=rowcount+1
    }
    names(pvallist)=rownames(currreal)
    enrichpvals[[count]]=pvallist
    count=count+1
  }
  enrichpvals
}

getNullCor=function(traitvec, RERmat, trimmedtree, genetrees, type="simperm", winR=NULL, winT=NULL){
  if(!type %in% c("simperm", "sim", "perm")){
    stop("type must be simperm, sim, or perm")
  }

  #get new vector
  if(type=="simperm"){
    vec=simpermvec(traitvec, trimmedtree)
  }
  if(type=="sim"){
    vec=simulatevec(traitvec, trimmedtree)
  }
  if(type=="perm"){
    vec=permutevec(traitvec)
  }

  #get enrich results
  paths=RERconverge::char2Paths(vec, genetrees)
  out=RERconverge::getAllCor(RERmat, paths, method="p", min.pos=0, winsorizeRER=winR, winsorizetrait = winT)
  out
}

getPermpval=function(permmatrix, realenrich){
  #returns vector of permpvals
  #sort real enrich
  realenrich=realenrich[order(rownames(realenrich)),]
  realpvals=realenrich$pval
  permpvals=c()
  count=1
  while(count<=nrow(realenrich)){
    numsig=sum(permmatrix[count,]<=realpvals[count])
    thisp=numsig/ncol(permmatrix)
    permpvals=c(permpvals, thisp)
    count=count+1
  }
  realenrich$permpval=permpvals
  realenrich
}

permutevec=function(namedvec){
  #returns permuted vec
  n=names(namedvec)
  vec=sample(namedvec)
  names(vec)=n
  vec
}

simulatevec=function(namedvec, treewithbranchlengths){
  #returns simulated vec
  #tree must be rooted and fully dichotomous
  #species in tree must match species in vec
  library("geiger")
  rm=ratematrix(treewithbranchlengths, namedvec)
  sims=sim.char(treewithbranchlengths, rm, nsim = 1)
  nam=rownames(sims)
  s=as.data.frame(sims)
  simulatedvec=s[,1]
  names(simulatedvec)=nam
  vec=simulatedvec
  vec
}

simpermvec=function(namedvec, treewithbranchlengths){
  #returns sim/perm vec
  #tree must be rooted and fully dichotomous
  #species in tree must match species in vec
  #simulate vector
  vec=simulatevec(namedvec, treewithbranchlengths)

  #assign real values to vec
  simsorted=sort(vec)
  realsorted=sort(namedvec)
  l=length(simsorted)
  c=1
  while(c<=l){
    simsorted[c]=realsorted[c]
    c=c+1
  }
  simsorted
}




