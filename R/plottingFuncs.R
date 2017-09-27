#comment everything out

require(dplyr)
require(ggplot2)

if(F){
multiplot = function(..., plotlist=NULL, file, cols=1, layout=NULL, widths=NULL, heights=NULL, flip=F) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  if(flip){
    layout=t(layout)
  }
  if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()



    if(!is.null(widths)){
      widths=widths
    }
    else{
      widths = unit(rep_len(1, ncol(layout)), "null")
    }
    if(!is.null(heights)){
      heights=heights
    }
    else{
      heights = unit(rep_len(1, nrow(layout)), "null")
    }
    show(widths)
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout), widths = widths, heights = heights )))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}


plotAllHists=function(resAll){
  plotlist=list()
  count=0;
  for(i in 1:length(resAll$res)){
    count=count+1
    p=resAll$res[[i]]$P
    s= sign(resAll$res[[i]]$Rho)
    ii=which(s!=0)
    p=p[ii]
    s=s[ii]
    plotlist[[count]]=plotHists(p, as.factor(s))+ggtitle(names(resAll$res)[i])
  }
  return(plotlist)

}

plotAllHistsWithPerm=function(resAll, resAllSim, resAllPerm=NULL, updown=T){
  plotlist=list()
  count=1;
  for(i in 1:length(resAll$res)){
    #  count=count+1
    realp=resAll$res[[i]]$P
    reals= sign(resAll$res[[i]][[1]])
    ii=which(reals!=0)
    realp=realp[ii]
    reals=reals[ii]


    simp=resAllSim[[i]]
    sims=sign(simp)
    simp=abs(simp)
    if(updown){
      show("UpDown")
      if(!is.null(resAllPerm)){
        permp=resAllPerm[[i]]
        perms=sign(permp)
        permp=abs(permp)

        vals=c(tmp1<-realp[reals>0], tmp2<-permp[perms>0], tmp3<-simp[sims>0])
        vals.grp=ordered(rep(c("Real", "Permuted","Simulated"), c(length(tmp1), length(tmp2),length(tmp3)), levels=c("Real", "Simulated","Permuted")))
      }
      else{
        vals=c(tmp1<-realp[reals>0], tmp3<-simp[sims>0])
        vals.grp=rep(c("Real","Simulated"), c(length(tmp1), length(tmp3)))
      }

      plotlist[[count]]=plotHists(vals, as.factor(vals.grp))+ggtitle(paste(names(resAll$res)[i], "Positive"))+theme_bw()+ theme(legend.position="none")
      count=count+1

      if(!is.null(resAllPerm)){
        vals=c(tmp1<-realp[reals<0], tmp2<-permp[perms<0], tmp3<-simp[sims<0])
        vals.grp=ordered(rep(c("Real", "Permuted","Simulated"), c(length(tmp1), length(tmp2),length(tmp3)), levels=c("Real", "Simulated","Permuted")))
      }
      else{
        vals=c(tmp1<-realp[reals<0], tmp3<-simp[sims<0])
        vals.grp=rep(c("Real","Simulated"), c(length(tmp1), length(tmp3)))
      }

      plotlist[[count]]=plotHists(vals, as.factor(vals.grp))+ggtitle(paste(names(resAll$res)[i], "Negative"))+theme_bw()
      count=count+1

    }

    else{
      show("All")
      if(!is.null(resAllPerm)){
        permp=resAllPerm[[i]]
        perms=sign(permp)
        permp=abs(permp)

        vals=c(tmp1<-realp[], tmp2<-permp[], tmp3<-simp[])
        vals.grp=ordered(rep(c("Real", "Permuted","Simulated"), c(length(tmp1), length(tmp2),length(tmp3)), levels=c("Real", "Simulated","Permuted")))
      }
      else{
        vals=c(tmp1<-realp[], tmp3<-simp[])
        vals.grp=rep(c("Real","Simulated"), c(length(tmp1), length(tmp3)))
      }
      show(count)
      plotlist[[count]]=plotHists(vals, as.factor(vals.grp))+ggtitle(paste(names(resAll$res)[i], "All"))+theme_bw()
      count=count+1
      show(count)
    }
  }
  return(plotlist)

}


plotHists=function(data, grp, ...){
  ii=which(is.na(data))
  if(length(ii)>0){
    data=data[-ii]
    grp=grp[-ii]
  }
  x=data.frame(data=data)
  x$grp=grp
  ggplot(x, aes(data, fill=grp))+geom_histogram(aes(y=..density..),position="dodge",  alpha=1,right=T, origin=0, ...)
}



plotContinuousChar=function(gene, treeObj, tip.vals, tip.vals.ref=NULL, rank=F, nlevels=8, type="c", col=NULL, residfun=residLO, useDiff=T){
  #get the tree projection
  tip.vals=tip.vals[!is.na(tip.vals)]
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  stopifnot(gene %in% names(treeObj$trees))
  tree=treeObj$trees[[gene]]
  stopifnot(!is.null(names(tip.vals)))
  both=intersect(tree$tip.label, names(tip.vals))

  stopifnot(length(both)>10)


  torm=setdiff(treeObj$masterTree$tip.label, both)
  tree=unroot(pruneTree(tree, both))
  allreport=treeObj$report[,both]
  ss=rowSums(allreport)
  iiboth=which(ss==length(both))


  ee=edgeIndexRelativeMaster(tree, treeObj$masterTree)
  ii= match(namePaths(ee,T), colnames(treeObj$paths))

  allbranch=treeObj$paths[iiboth,ii]

  allbranch=scaleMat_c(allbranch)
  show(sum(is.na(allbranch)))
  nv=projection(t(allbranch), method="AVE", returnNV = T)

  proj=residfun(tree$edge.length, nv)
  treeChar=edgeVars(tree, tip.vals, useDiff=useDiff)


  nn=nameEdges(tree)
  nn[nn!=""]=speciesNames[nn[nn!=""], ]
  par(mfrow=c(2,2), mai=rep(0.7,4))

  plotWtext(sqrt(nv), sqrt(tree$edge.length), xlab="char", ylab="Gene branch length", labels = nn)
  if(!is.null(tip.vals.ref)){
    treeCharRef=edgeVars(tree, tip.vals.ref, useDiff=useDiff)
    proj=resid(rbind(proj), model.matrix(~1+treeCharRef$edge.length))[1,]
  }
  plotWtext(treeChar$edge.length, proj, xlab="char", ylab="Gene branch length", labels = nn)
  stat=cor.test(treeChar$edge.length, proj, method="s")
  mtext(paste0("r=", round(stat$estimate,2), ";  p-value=", format.pval(stat$p.value)), side = 3, line=1)



  tree$tip.label=speciesNames[tree$tip,1]
  col=colorpanel(nlevels, "blue", "red", "yellow3")


  vals=treeChar$edge.length
  if(rank){
    vals=rank(vals)
  }
  plot.phylo(tree, use.edge.length = F,type=type,edge.color=col[cut(vals, nlevels)], edge.width=3.5, lab4ut="axial", cex=0.6)
  title("character")

  vals=proj
  if(rank){
    vals=rank(vals)
  }
  plot.phylo(tree, use.edge.length = F,type=type,edge.color=col[cut(vals, nlevels)], edge.width=3.5, lab4ut="axial", cex=0.6)
  title("projection")



}

treePlot=function(tree, vals=NULL,rank=F, nlevels=8, type="c", col=NULL){
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  if(is.null(vals)){
    vals=tree$edge.length
  }
  vals=as.numeric(vals)
  if(rank){
    vals=rank(vals)
  }
  layout(matrix(c(1,2), ncol=1),heights=c(10,2))
  if(is.null(col)){
    col=colorpanel(nlevels, "blue", "red")
  }
  tree$tip.label=speciesNames[tree$tip,1]
  plot.phylo(tree, use.edge.length = F,type=type,edge.color=col[cut(vals, nlevels)], edge.width=3.5, lab4ut="axial", cex=0.6)


  min.raw <- min(vals, na.rm = TRUE)
  max.raw <- max(vals, na.rm = TRUE)
  z <- seq(min.raw, max.raw, length = length(col))

  par(mai=c(1,0.5,0,0.5))
  image(z = matrix(z, ncol = 1), col = col, breaks = seq(min.raw, max.raw, length.out=nlevels+1),
        xaxt = "n", yaxt = "n")
  #par(usr = c(0, 1, 0, 1))
  lv <- pretty(seq(min.raw, max.raw, length.out=nlevels+1))
  scale01 <- function(x, low = min(x), high = max(x)) {
    x <- (x - low)/(high - low)
    x
  }
  xv <- scale01(as.numeric(lv), min.raw, max.raw)
  axis(1, at = xv, labels = lv)
}

}

plotRers <- function(rermat, index, phenv, method = 's'){
     e1 = rermat[index,][!is.na(rermat[index,])]     
     colids = !is.na(rermat[index,])
     phenvid = phenv[colids]          
     e1plot <- e1
     if(exists('speciesNames')){
          names(e1plot) <- speciesNames[names(e1),]
     }
     fgdcor = getAllCor(rermat[index,,drop=F],phenv, method = method)
     if(is.numeric(index)){
         gen = rownames(rermat)[index]  
     }else{
         gen = index
     }
     plottitle = paste0(gen, ': rho = ',round(fgdcor$Rho,4),', p = ',round(fgdcor$P,4))
     print(plottitle)
     names(e1plot)[is.na(names(e1plot))]=""
     fgd = setdiff(names(e1plot)[phenvid == 1],"")
     df <- data.frame(species = names(e1plot), rer = e1plot, stringsAsFactors=FALSE) %>%
          mutate(mole = as.factor(ifelse(names(e1plot) %in% fgd,2,1)))
     ll=c(min(df$rer)*1.1, max(df$rer)+0.2)     
     g  <- ggplot(df, aes(x = rer, y=factor(species), col=mole, label=species)) + scale_size_manual(values=c(3,3))+ geom_point(aes(size=mole))+
          scale_color_manual(values = c("deepskyblue3", "brown1"))+
          scale_x_continuous(limits=ll)+
          #scale_x_continuous(expand = c(.1,.1))+
          geom_text(hjust=1, size=5)+
          ylab("Branches")+
          xlab("relative rate")+
          ggtitle(plottitle)+
          geom_vline(xintercept=0, linetype="dotted")+
          theme(axis.ticks.y=element_blank(),axis.text.y=element_blank(),legend.position="none",
                panel.background = element_blank(),
                axis.text=element_text(size=18,face='bold',colour = 'black'),
                axis.title=element_text(size=24,face="bold"),
                plot.title= element_text(size = 24, face = "bold"))+
          theme(axis.line = element_line(colour = 'black',size = 1))+
                theme(axis.line.y = element_blank())
          
     print(g)
}