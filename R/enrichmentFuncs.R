#comment everything out for now
if(T){
wilcoxGMTall=function(vals, annotList){
  reslist=list()
  for ( n in names(annotList)){
    reslist[[n]]=wilcoxGMT(vals, annotList[[n]])
    message(paste0(nrow(reslist[[n]]), " results for annotation set ", n))
  }
  reslist
}

wilcoxGMT=function(vals, gmt, simple=F, use.all=F, num.g=10,genes=NULL){
     #vals = resstat
     #gmt = annotlist
     #simple = F
     #use.all = F
     #num.g = 10
     #genes = NULL
  vals=vals[!is.na(vals)]
  if(is.null(genes)){
    genes=unique(unlist(gmt$genesets))
  }
  out=matrix(nrow=length(gmt$genesets), ncol=3)
  rownames(out)=gmt$geneset.names
  colnames(out)=c("stat", "pval", "genes")
  out=as.data.frame(out)
  genes=intersect(genes, names(vals))
  show(length(intersect(genes, names(vals))))
  show(length(genes))
  #vals=rank(vals[genes])
  for( i in 1:nrow(out)){

    curgenes=intersect(genes,gmt$genesets[[i]])

    bkgenes=setdiff(genes, curgenes)

    if (length(bkgenes)==0 || use.all){
      bkgenes=setdiff(names(vals), curgenes)
    }
    if(length(curgenes)>=num.g & length(bkgenes)>2){
      if(!simple){
        res=wilcox.test(x = vals[curgenes], y=vals[bkgenes])

        out[i, 1:2]=c(res$statistic/(length(bkgenes)*length(curgenes)), res$p.value)
      }
      else{
        valsr=rank(vals[genes])
        out[i, 1:2]=simpleAUCgenesRanks(valsr[curgenes],valsr[bkgenes])

      }
      if (out[i,1]>0.5){
        oo=order(vals[curgenes], decreasing = T)
        granks=rank(-vals)
      }
      else{
        oo=order(vals[curgenes], decreasing = F)
        granks=rank(vals)
      }
      # show(vals[curgenes][oo])

      nn=paste(curgenes[oo],round((granks[curgenes])[oo],2),sep=':' )
      out[i,3]=length(curgenes)
      out[i,4]=paste(nn, collapse = ", ")

    }

  }
  # hist(out[,2])
  out[,1]=out[,1]-0.5
  #out=cbind(out,BH(out[,2]))
  #out=cbind(out,p.adjust(out[,2], method = 'BH'))
  #colnames(out)[ncol(out)]="p.adj"

  # out=out[, c(1,2,3,5,4)]
  out=out[!is.na(out[,2]),]
  #out=out[out[,3]<0.2,]
  out=out[order(-abs(out[,1])),]
}


simpleAUCmat=function(lab, value){
  value=t(apply(rbind(value),1,rank))
  posn=sum(lab>0)
  negn=sum(lab<=0)
  if(posn<2||negn<2){
    auc=rep(NA, nrow(value))
    pp=rep(NA, nrow(value))
  }
  else{
    stat=apply(value[,lab>0, drop=F],1,sum)-posn*(posn+1)/2

    auc=stat/(posn*negn)
    mu=posn*negn/2
    sd=sqrt((posn*negn*(posn+negn+1))/12)
    stattest=apply(cbind(stat, posn*negn-stat),1,max)
    pp=(2*pnorm(stattest, mu, sd, lower.tail = F))
  }
  return(list(auc=auc, pp=pp))
}


myAUC<-function(labels, values){
  ii=which(!is.na(values))
  values=values[ii]
  labels=labels[ii]
  posii=which(labels>0)
  negii=which(labels<=0)
  posn=length(posii)
  negn=length(negii)
  posval=values[posii]
  negval=values[negii]
  res=wilcox.test(posval, negval, alternative="greater", conf.int=TRUE, exact=F, correct=F);

  myres=list()
  myres$low=res$conf.int[1]
  myres$high=res$conf.int[2]

  myres$auc=(res$statistic)/(posn*negn)
  myres$pval=res$p.value
  return(myres)
}

simpleAUCgenesRanks=function(pos, neg){

  posn=length(pos)
  negn=length(neg)
  stat=sum(pos)-posn*(posn+1)/2
  auc=stat/(posn*negn)
  mu=posn*negn/2
  sd=sqrt((posn*negn*(posn+negn+1))/12)
  stattest=apply(cbind(stat, posn*negn-stat),1,max)
  pp=(2*pnorm(stattest, mu, sd, lower.tail = F))

  return(c(auc,pp))
}



fastwilcoxGMTall=function(vals, annotList, ...){
  reslist=list()
  for ( n in names(annotList)){
    reslist[[n]]=fastwilcoxGMT(vals, annotList[[n]], ...)
    message(paste0(nrow(reslist[[n]]), " results for annotation set ", n))
  }
  reslist
}

fastwilcoxGMT=function(vals, gmt, simple=T, use.all=F, num.g=10,genes=NULL, outputGeneVals=F, order=F){
  vals=vals[!is.na(vals)]
  if(is.null(genes)){
    genes=unique(unlist(gmt$genesets))
  }
  out=matrix(nrow=length(gmt$genesets), ncol=5)
  rownames(out)=gmt$geneset.names
  colnames(out)=c("stat", "pval", "p.adj","num.genes", "gene.vals")
  out=as.data.frame(out)
  genes=intersect(genes, names(vals))

  valsr=rank(vals[genes])
  numg=length(vals)+1
  valsallr=rank(vals)
  for( i in 1:nrow(out)){

    curgenes=intersect(genes,gmt$genesets[[i]])

    bkgenes=setdiff(genes, curgenes)

    if (length(bkgenes)==0 || use.all){
      bkgenes=setdiff(names(vals), curgenes)
    }
    if(length(curgenes)>=num.g & length(bkgenes)>2){
      if(!simple){
        res=wilcox.test(x = vals[curgenes], y=vals[bkgenes], exact=F)

        out[i, 1:2]=c(res$statistic/(length(bkgenes)*length(curgenes)), res$p.value)
      }
      else{

        out[i, 1:2]=simpleAUCgenesRanks(valsr[curgenes],valsr[bkgenes])

      }
      out[i,"num.genes"]=length(curgenes)
      if(outputGeneVals){
        if (out[i,1]>0.5){
          oo=order(vals[curgenes], decreasing = T)
          granks=numg-valsallr[curgenes]
        }
        else{
          oo=order(vals[curgenes], decreasing = F)
          granks=valsallr[curgenes]
        }


        nn=paste(curgenes[oo],round((granks[curgenes])[oo],2),sep=':' )
        out[i,"gene.vals"]=paste(nn, collapse = ", ")
      }
    }

  }
  # hist(out[,2])
  out[,1]=out[,1]-0.5
  out[, "p.adj"]=p.adjust(out[,2], method="BH")

  out=out[!is.na(out[,2]),]
  if(order){
    out=out[order(-abs(out[,1])),]
  }
  out
}

}

#adopted from GSA package, gsa.read.gmt()
read.gmt <- function (filename) 
{
    a = scan(filename, what = list("", ""), sep = "\t", quote = NULL, 
        fill = T, flush = T, multi.line = F)
    geneset.names = a[1][[1]]
    geneset.descriptions = a[2][[1]]
    dd = scan(filename, what = "", sep = "\t", quote = NULL)
    nn = length(geneset.names)
    n = length(dd)
    ox = rep(NA, nn)
    ii = 1
    for (i in 1:nn) {
        cat(i)
        while ((dd[ii] != geneset.names[i]) | (dd[ii + 1] != 
            geneset.descriptions[i])) {
            ii = ii + 1
        }
        ox[i] = ii
        ii = ii + 1
    }
    genesets = vector("list", nn)
    for (i in 1:(nn - 1)) {
        cat(i, fill = T)
        i1 = ox[i] + 2
        i2 = ox[i + 1] - 1
        geneset.descriptions[i] = dd[ox[i] + 1]
        genesets[[i]] = dd[i1:i2]
    }
    geneset.descriptions[nn] = dd[ox[nn] + 1]
    genesets[[nn]] = dd[(ox[nn] + 2):n]
    out = list(genesets = genesets, geneset.names = geneset.names, 
        geneset.descriptions = geneset.descriptions)
    class(out) = "GSA.genesets"
    return(out)
}
