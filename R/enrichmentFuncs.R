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
}
