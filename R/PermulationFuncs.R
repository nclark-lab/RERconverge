
#'Combines batches of permulations
#' @param permdat1 Batch of permulations output from \code{\link{getPermsContinuous}}
#' @param permdat2 Batch of permulations output from \code{\link{getPermsContinuous}}
#' @param enrich Default T. Specifies if both `permdat1` and `permdat2` contain enrichment permulations
#' @return Combined permulations
#' @export
combinePermData=function(permdat1, permdat2, enrich=T){
  #combine results
  allcorP=merge(permdat1$corP,permdat2$corP,by="row.names")
  allcorRho=merge(permdat1$corRho,permdat2$corRho,by="row.names")
  allcorStat=merge(permdat1$corStat,permdat2$corStat,by="row.names")
  rownames(allcorP)=allcorP[,1]
  allcorP[,1]=NULL
  rownames(allcorRho)=allcorRho[,1]
  allcorRho[,1]=NULL
  rownames(allcorStat)=allcorStat[,1]
  allcorStat[,1]=NULL

  if(enrich){
    #combine enrichment
    allenrichP=permdat1$enrichP
    allenrichP2=permdat2$enrichP
    allenrichStat=permdat1$enrichStat
    allenrichStat2=permdat2$enrichStat
    g=1
    while(g<=length(allenrichP)){
      allenrichP[[g]]=merge(allenrichP[[g]], allenrichP2[[g]],by="row.names")
      rownames(allenrichP[[g]])=allenrichP[[g]][,1]
      allenrichP[[g]][,1]=NULL

      allenrichStat[[g]]=merge(allenrichStat[[g]], allenrichStat2[[g]],by="row.names")
      rownames(allenrichStat[[g]])=allenrichStat[[g]][,1]
      allenrichStat[[g]][,1]=NULL

      g=g+1
    }
  }

  if(enrich){
    data=vector("list", 5)
    data[[1]]=allcorP
    data[[2]]=allcorRho
    data[[3]]=allcorStat
    data[[4]]=allenrichP
    data[[5]]=allenrichStat
    names(data)=c("corP", "corRho", "corStat", "enrichP", "enrichStat")
  }else{
    data=vector("list", 3)
    data[[1]]=allcorP
    data[[2]]=allcorRho
    data[[3]]=allcorStat
    names(data)=c("corP", "corRho", "corStat")
  }
  data
}

#'Calculates permuted correlation and enrichment statistics
#' @param numperms An integer number of permulations
#' @param traitvec A named phenotype vector
#' @param RERmat An RER matrix calculated using \code{\link{getAllResiduals}}
#' @param annotlist Pathway annotations
#' @param trees treesObj from \code{\link{readTrees}}
#' @param mastertree A rooted, fully dichotomous tree derived from the treesObj master tree from \code{\link{readTrees}}.  Must not contain species not in traitvec
#' @param calculateenrich A boolean variable indicating if null permulation p-values for enrichment statistics
#' @param type One of "simperm", "sim", or "perm" for permulations, simulations, or permutations, respectively
#' @param winR Integer winzorization value for RERs
#' @param winT Integer winzorization value for trait
#' @param method statistical method to use for correlations
#' @param min.pos minimum foreground species - should be set to 0
#' @note  winsorize is in terms of number of observations at each end, NOT quantiles
#' @return A list object with enrichment statistics, correlation p-val, rho, and correlation effect size
#' @export
getPermsContinuous=function(numperms, traitvec, RERmat, annotlist, trees, mastertree, calculateenrich=T, type="simperm", winR=3, winT=3, method="p", min.pos=0){

  #get real enrich and cors
  realpaths=RERconverge::char2Paths(traitvec, trees)
  realresults=RERconverge::getAllCor(RERmat, realpaths, method = method, min.pos = min.pos, winsorizeRER = winR, winsorizetrait=winT)
  realstat=sign(realresults$Rho)*-log10(realresults$P)

  #make enrich list/matrices to fill
  permPvals=data.frame(matrix(ncol=numperms, nrow=nrow(realresults)))
  rownames(permPvals)=rownames(realresults)
  permRhovals=data.frame(matrix(ncol=numperms, nrow=nrow(realresults)))
  rownames(permRhovals)=rownames(realresults)
  permStatvals=data.frame(matrix(ncol=numperms, nrow=length(realstat)))
  rownames(permStatvals)=rownames(realresults)


  if(calculateenrich){
    realenrich=fastwilcoxGMTall(na.omit(realstat), annotlist, outputGeneVals=F)

    #sort real enrichments
    groups=length(realenrich)
    c=1
    while(c<=groups){
      current=realenrich[[c]]
      realenrich[[c]]=current[order(rownames(current)),]
      c=c+1
    }
    #make matrices to fill
    permenrichP=vector("list", length(realenrich))
    permenrichStat=vector("list", length(realenrich))
    c=1
    while(c<=length(realenrich)){
      newdf=data.frame(matrix(ncol=numperms, nrow=nrow(realenrich[[c]])))
      rownames(newdf)=rownames(realenrich[[c]])
      permenrichP[[c]]=newdf
      permenrichStat[[c]]=newdf
      c=c+1
    }
  }


  counter=1;
  while(counter<=numperms){

    print(paste("running permutation: ", counter))

    #get correlation results
    out=getNullCor(traitvec, RERmat, mastertree, trees, type = type, winR=winR, winT=winT)
    stat=sign(out$Rho)*-log10(out$P)

    permPvals[,counter]=out$P
    permRhovals[,counter]=out$Rho
    permStatvals[,counter]=stat

    if(calculateenrich){
      enrich=fastwilcoxGMTall(na.omit(stat), annotlist, outputGeneVals=F)
      #sort and store enrichment results
      groups=length(enrich)
      c=1
      while(c<=groups){
        current=enrich[[c]]
        enrich[[c]]=current[order(rownames(current)),]
        permenrichP[[c]][,counter]=enrich[[c]]$pval
        permenrichStat[[c]][,counter]=enrich[[c]]$stat
        c=c+1
      }
    }
    counter=counter+1
  }

  if(calculateenrich){
    data=vector("list", 5)
    data[[1]]=permPvals
    data[[2]]=permRhovals
    data[[3]]=permStatvals
    data[[4]]=permenrichP
    data[[5]]=permenrichStat
    names(data)=c("corP", "corRho", "corStat", "enrichP", "enrichStat")
  }else{
    data=vector("list", 3)
    data[[1]]=permPvals
    data[[2]]=permRhovals
    data[[3]]=permStatvals
    names(data)=c("corP", "corRho", "corStat")
  }
  data
}

#'Performs enrichment statistic permulations using existing gene correlation permulations
#' @param corperms Gene correlation permulations from \code{\link{getPermsContinuous}}
#' @param realenrich Pathway enrichment results using observed phenotype obtained from `correlateWithContinuousPhenotype`
#' @param annotlist Pathway annotations
#' @return Full null permulation statistics and p-values for gene correlations and pathway enrichment
#' @export
getEnrichPermsContinuous=function(corperms, realenrich, annotlist){
  numperms=ncol(corperms$corP)
  #sort real enrichments
  groups=length(realenrich)
  c=1
  while(c<=groups){
    current=realenrich[[c]]
    realenrich[[c]]=current[order(rownames(current)),]
    c=c+1
  }
  #make matrices to fill
  permenrichP=vector("list", length(realenrich))
  permenrichStat=vector("list", length(realenrich))
  c=1
  while(c<=length(realenrich)){
    newdf=data.frame(matrix(ncol=numperms, nrow=nrow(realenrich[[c]])))
    rownames(newdf)=rownames(realenrich[[c]])
    permenrichP[[c]]=newdf
    permenrichStat[[c]]=newdf
    c=c+1
  }

  count=1
  while(count<=numperms){
    print(paste("running permutation: ", count))
    statvec=setNames(corperms$corStat[,count], rownames(corperms$corP))
    statvec=na.omit(statvec)
    enrich=fastwilcoxGMTall(statvec, annotlist, outputGeneVals=F)
    #sort and store enrichment results
    groups=length(enrich)
    c=1
    while(c<=groups){
      current=enrich[[c]]
      enrich[[c]]=current[order(rownames(current)),]
      permenrichP[[c]][,count]=enrich[[c]]$pval
      permenrichStat[[c]][,count]=enrich[[c]]$stat
      c=c+1
    }
    count=count+1
  }
  data=vector("list", 5)
  data[[1]]=corperms$corP
  data[[2]]=corperms$corRho
  data[[3]]=corperms$corStat
  data[[4]]=permenrichP
  data[[5]]=permenrichStat
  names(data)=c("corP", "corRho", "corStat", "enrichP", "enrichStat")
  data
}

#'Calculates enrichment permutation pvals from output of \code{\link{getPermsContinuous}}
#' @param realenrich Real enrichment statistics from \code{\link{fastwilcoxGMTall}}
#' @param permvals output from \code{\link{getPermsContinuous}}
#' @return A list object with vectors of permulation p-values
#' @export
permpvalenrich=function(realenrich, permvals){

  #sort real enrichments
  groups=length(realenrich)
  c=1
  while(c<=groups){
    current=realenrich[[c]]
    realenrich[[c]]=current[match(rownames(permvals$enrichStat[[c]]), rownames(current)),]
    c=c+1
  }

  permenrich=permvals$enrichStat
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
      if(is.na(currreal[rowcount,]$stat)){
        pval=lessnum/denom
      }else{
        lessnum=sum(abs(currenrich[rowcount,])>abs(currreal[rowcount,]$stat), na.rm=T)
        denom=sum(!is.na(currenrich[rowcount,]))
        pval=lessnum/denom
      }
      pvallist=c(pvallist, pval)
      rowcount=rowcount+1
    }
    names(pvallist)=rownames(currreal)
    enrichpvals[[count]]=pvallist
    count=count+1
  }
  enrichpvals
}

#'Calculates correlation permutation pvals from output of \code{\link{getPermsContinuous}}
#' @param realcor Real enrichment statistics from \code{\link{fastwilcoxGMTall}}
#' @param permvals output from \code{\link{getPermsContinuous}}
#' @return A vector with permulation p-values
#' @export
permpvalcor=function(realcor, permvals){

  permcor=permvals$corRho
  realstat=realcor$Rho
  names(realstat)=rownames(realcor)

  permcor=permcor[match(names(realstat), rownames(permcor)),]

  permpvals=vector(length=length(realstat))
  names(permpvals)=names(realstat)
  count=1
  while(count<=length(realstat)){
    if(is.na(realstat[count])){
      permpvals[count]=NA
    }else{
      num=sum(abs(permcor[count,])>abs(realstat[count]), na.rm=T)
      denom=sum(!is.na(permcor[count,]))
      permpvals[count]=num/denom
    }
    count=count+1
  }
  permpvals
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
  paths=RERconverge::char2Paths(na.omit(vec), genetrees)
  out=RERconverge::getAllCor(RERmat, paths, method="p", min.pos=0, winsorizeRER=winR, winsorizetrait = winT)
  out
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


