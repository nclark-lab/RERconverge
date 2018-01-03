
#' @keywords  internal
projection <- function(protein, rna=rvector, method=c("RNA","AVE","PCA"), returnNV=F)
{


  ###projection(ribosomal RNA)###
  if (match.arg(method) == "RNA") {
    rvector <- rna[,1]
    rvector.sc <- rvector/sqrt(sum(rvector^2))
    normv=rvector.sc
    result <- as.matrix(protein) - rvector.sc %*% t(rvector.sc) %*% as.matrix(protein)
  }

  ###projection(average vector)###
  if (match.arg(method) == "AVE") {
    n <- ncol(protein)
    m <- nrow(protein)
    protein.sc <- matrix(0,m,n)
    for (i in 1:n) {
      unit <- protein[,i] / sqrt(sum(protein[,i]^2,na.rm=TRUE))
      protein.sc[,i] <- unit
    }
    mean.clean = function(x){mm=mean(x,na.rm=TRUE); mm }
    #  average <- apply(protein.sc, 1, mean.clean)
    average=rowMeans(protein.sc, na.rm=T)
    avector.sc <- average/sqrt(sum(average^2))
    #    result <- average/sqrt(sum(average^2))
    normv=avector.sc
    Id=diag(nrow(protein))
    #result=(Id-avector.sc %*% t(avector.sc))%*%protein
    result <- as.matrix(protein) - avector.sc %*% t(avector.sc) %*% as.matrix(protein)
  }

  ###projection(PCA)###
  if (match.arg(method) == "PCA") {
    pvector <- as.vector(prcomp(protein)$x[,1])
    print(pvector)
    pvector.sc <- pvector/sqrt(sum(pvector^2))
    normv=pvector.sc
    result <- as.matrix(protein) - pvector.sc %*% t(pvector.sc) %*% as.matrix(protein)
  }
  if(returnNV){
    normv
  }
  else{
    result
  }
}
#' @keywords internal
resid=function(dat, lab){
  if(is.null(dim(dat))){
    dat=rbind(dat)
  }
  if (is.null(dim(lab))){
    mod=model.matrix(~1+lab);
  }
  else{
    mod=lab
  }
  n=dim(dat)[2]
  Id=diag(n)
  resid=dat %*% (Id - mod %*% solve(t(mod) %*% mod) %*%
                   t(mod))
}
#' @keywords  internal
naresid=function(data, X, useZero=F, weights=NULL, covar=NULL){
  if(is.vector(X)){
    mod=model.matrix(~1+as.vector(X));
  }
  else{
    mod=X
  }
  #show(mod)
  resid=matrix(nrow=nrow(data), ncol=ncol(data))
  if(useZero){
    resid[]=0;
  }
  else{
    resid[]=NA;
  }
  for ( i in 1:nrow(data)){

    ii=which(!is.na(data[i,]))
    if(length(ii)>2){
      iiinv=which(is.na(data[i,]))

      dat=data[i,ii,drop=F]

      modtmp=mod[ii,]
      if(!is.null(weights)){
        W=diag(weights[i,ii])
      }
      else{
        W=NULL
      }
      if (!is.null(covar)){
        if(!is.null(W)){
          W=W%*%covar
        }
        else{
          W=covar
        }
      }
      n=dim(dat)[2]
      Id=diag(n)
      if(!is.null(W)){

        coeff=dat%*%W%*%modtmp %*% solve(t(modtmp) %*% W %*% modtmp)
        #check there is no error with lm
        # lmres=lm(t(dat)~0+modtmp, weights = diag(W))
        #  show(coeff)
        #  show(coefficients(lmres))
        resid[i, ii] = dat -(coeff %*% t(modtmp))
        resid[i, ii]=resid[i,ii]*sqrt(diag(W))


      }
      else{
        resid[i,ii]=dat %*% (Id - modtmp %*% solve(t(modtmp) %*% modtmp) %*%
                               t(modtmp))
      }



    }
    else{
      # message("Cannot compute residuals")

    }

  }
  rownames(resid)=rownames(data)
  colnames(resid)=colnames(data)
  resid
}




#' @keywords  internal
projectRLM=function(data, nv){
  resid=matrix(nrow=nrow(data), ncol=ncol(data))
  rownames(resid)=rownames(data)
  colnames(resid)=colnames(data)
  for ( i in 1:nrow(data)){
    print(i)
    ii=which(!is.na(data[i,]))
    iir=which(is.na(data[i,]))
    y=data[i,ii]
    nvtmp=as.vector(nv[ii])
    lmres=rlm(y~1+nvtmp, psi=psi.huber)
    resid[i,ii]=lmres$residuals
    resid[i, iir]=0
  }
  resid
}
