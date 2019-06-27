#' RERconverge
#'
#'
#' @docType package
#' @author
#' @import Rcpp
#' @importFrom Rcpp evalCpp
#' @useDynLib RERconverge
#' @name RERconverge
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

treePlot=function(tree, vals=NULL,rank=F, nlevels=5, type="c", col=NULL){
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
treePlot=function(tree, vals=NULL,rank=F, nlevels=5, type="c", col=NULL){
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
  #tree$tip.label=speciesNames[tree$tip,1]
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

#Version of tree plotting function with more options for colors and tip labels
#Used for plotting trees for PON1 manuscript
#Issue with invalid graphics state after running: Solve by only plotting legend if
#margins are not too large (how to test?)

#' Plot `tree` with branch labels colored in a heatmap based on values in `vals`
#'
#' @param tree. A phylo object, used for topology
#' @param vals. Values to use for heatmap to color branches, in same order as edges of tree
#' @param rank. Whether to plot ranks instead of values
#' @param nlevels. How many colors to use in the heatmap
#' @param type. Type for plot.phylo
#' @param col. Vector of user-defined colors
#' @param maintitle. Main title label
#' @param useedge. Whether to use edge lengths in `tree` for plotting
#' @param doreroot. Whether to re-root the tree before  plotting
#' @param rerootby. If re-rooting, what to use to root the tree (verify by checking against unrooted plot)
#' @param useSpecies. A vector of species to include in the plot (verify by checking against full set)
#' @param species.names. Data fram for converting names in `tree` to names to be plotted
#' (row names should be tip labels of tree, and first column should contain the
#' corresponding desired tip labels)
#' @param speclist1. A vector of tip labels to highlight in bold, blue text
#' @param speclist2. A vector of tip labels to which to add an asterisk
#' @param aligntip. Whether to align tip labels (default FALSE)
#' @param colpan1. Color for lowest value in heatmap (default blue)
#' @param colpan2. Color for highest value in heatmap (default rose)
#' @param colpanmid. Color for middle value in heatmap (default gray)
#' @param plotspecies. A vector of tip labels to display on the tree (the remainder will be masked,
#' but the corresponding tips will be plotted on the tree)
#' @param edgetype. Vector of line type for edges.
#' @param textsize. cex value for tip labels (default 0.6)
#' @param colbarlab. Label for the color bar legend
#' @param splist2sym. A value within species names to display as a symbol
#' @param dolegend. Whether to display the heatmap legend.
#' @param nacol. Color to display for any edges with length NA
#' @param figwid. Adjust x limits of plot.phylo by 1/figwid. May be related to figure width
#' (requires some optimization)
#' @param .... further arguments to be passed to `plot` or to `plot.phylo`
#' @return Plots a cladogram of `tree` with branch colors determined by `vals`
#' @export

treePlotNew=function(tree, vals=NULL, rank=F, nlevels=9, type="c", col=NULL,
                     maintitle= NULL, useedge=F, doreroot=F, rerootby=NULL, useSpecies=NULL,
                     species.names=NULL, speclist1=NULL, speclist2=NULL, aligntip=F,
                    colpan1=rgb(0,119,187,maxColorValue=255),
                    colpan2=rgb(204,51,17,maxColorValue=255),
                    colpanmid=rgb(187,187,187,maxColorValue=255),plotspecies=NULL,
                    edgetype=NULL,textsize=0.6,colbarlab="",splist2sym="psi",
                    dolegend=T,nacol=rgb(0,0,0),figwid=10,...){
 #bold speclist1, star speclist2
 #reroot before plotting to match up vals
  require(gplots)
  if(is.null(vals)){
    vals=tree$edge.length
  }
  vals=as.numeric(vals)
  if (rank) {
    vals = rank(vals)
  }
  faketree=tree
  faketree$edge.length=vals
  #layout(matrix(c(1,2), ncol=1),heights=c(10,2))
  #layout(matrix(c(1,2), ncol=1),heights=c(7,1))
  if(is.null(col)){
    if (is.null(colpanmid)) {
      col=colorpanel(nlevels, colpan1, colpan2)
    } else {
      col=colorpanel(nlevels, colpan1, colpanmid, colpan2)
    }
  }
  if (!is.null(useSpecies)) {
    message("Pruning may influence color plotting; greater accuracy with branches displayed as NA")
    tree=pruneTree(tree,useSpecies)
    faketree=pruneTree(tree,useSpecies)
    vals=faketree$edge.length
  }
  if (doreroot) {
    message("Rooting may influence color plotting; greater accuracy with original topology")
    rerootby=intersect(rerootby,tree$tip.label)
    if (length(rerootby) > 0) {
        tree=root(tree,rerootby,resolve.root=T)
        #want to distribute vals along branches instead
        faketree=root(faketree,rerootby,resolve.root=T)
        vals=faketree$edge.length
    } else {
      print("No species in rerootby in tree! Leaving tree unrooted.")
    }
  }
  if (!is.null(species.names)) {
    for(s in 1:length(tree$tip.label)) {
      if (tree$tip.label[s] %in% row.names(species.names)) {
        tree$tip.label[s] <- species.names[,1][which(row.names(species.names)==tree$tip.label[s])]
        }
    }
  }
  pfonts <- c(rep(1,length(tree$tip.label)))
  tipcol <- c(rep("black",length(tree$tip.label)))
  if (!is.null(speclist1)) {
    pfonts[which(tree$tip.label %in% speclist1)] <- 2
    tipcol[which(tree$tip.label %in% speclist1)] <- "blue"
  }
  if (! is.null(plotspecies)) { #only plot certain species names by making others white
    tipcol[which(tree$tip.label %in% plotspecies == FALSE)] <- "white"
  }
  if (!is.null(speclist2)) {
    toadd <- which(tree$tip.label %in% speclist2)
    tree$tip.label[toadd] <- str_replace_all(tree$tip.label[toadd]," ","~")
    tree$tip.label[toadd] <- as.expression(parse(text=paste(tree$tip.label[toadd],"~",splist2sym,sep="")))
  }
  if (is.null(edgetype)) {
    edgetype = c(rep(1,length(vals))) #does not play well with rerootby
  }
  calcoff <- quantile(tree$edge.length[tree$edge.length>0],0.25,na.rm=T)
  if (min(vals,na.rm=T)>=0) {
    forbreaks = quantile(vals, probs = seq(0,1, length.out = nlevels+1), na.rm=T)
  } else {
    #separately estimate quantiles below and above 0
    negvals = vals[vals < 0]
    posvals = vals[vals >= 0]
    negbreaks = quantile(negvals, probs = seq(0,1, length.out = nlevels+1),
                         na.rm=T)[seq(1,nlevels+1,2)]
    posbreaks = quantile(posvals, probs = seq(0,1, length.out = nlevels+1),
                         na.rm=T)[seq(nlevels %% 2 + 1,nlevels+1,2)]
    forbreaks = unique(c(negbreaks,posbreaks))
  }
  eccalc = col[cut(vals, breaks = forbreaks, include.lowest = T, right = T)]
  edgetype[which(is.na(eccalc))] = 3 #set na branches to dashed
  eccalc[which(is.na(eccalc))] = nacol #set na branches to some other color (black?)

  #Plotting horizontally messes up the margins, so re-set x limits *before* setting layout
  #dev.new(width=10,height=9)
  pdf()
  dev.control('enable')
  layout(matrix(c(1,2), nrow=1),widths=c(5,1)) #switching to horizontal
  par(mar=c(2,1,0,0.2)+0.1)
  par(omi=c(0,0,0,0.0001))
  forx = plot.phylo(tree, use.edge.length = useedge,type=type,
                       edge.width=4, edge.lty=edgetype,lab4ut="axial", cex=textsize,
                       align.tip.label=aligntip,font=pfonts,label.offset=calcoff,
                       no.margin=T,plot=F,...)
  #print(forx$x.lim)
  dev.off()
  #dev.new(width=10,height=9)
  layout(matrix(c(1,2), nrow=1),widths=c(5,1)) #switching to horizontal
  par(mar=c(2,1,0,0.2)+0.1)
  #oldpar = par()
  #par(omi=c(0,0,0,0))
  par(omi=c(0,0,0,0.0001))
  #plotobj = plot.phylo(tree, use.edge.length = useedge,type=type,edge.color=col[cut(vals, nlevels)], edge.width=4, edge.lty=edgetype,lab4ut="axial", cex=textsize, align.tip.label=aligntip,font=pfonts,label.offset=calcoff,
  #tip.color=tipcol,no.margin=T,plot=T, main = maintitle)

  #Plotting horizontally messes up the margins, so re-set x limits
  plotobj = plot.phylo(tree, use.edge.length = useedge,type=type,
                       edge.color=eccalc,
                       edge.width=4, edge.lty=edgetype,lab4ut="axial", cex=textsize,
                       align.tip.label=aligntip,font=pfonts,label.offset=calcoff,
                       tip.color=tipcol,no.margin=T,plot=T, x.lim=forx$x.lim/figwid,main = maintitle,...)
  #par(mar=oldpar)
  if (dolegend) {
    min.raw <- min(vals, na.rm = TRUE)
    max.raw <- max(vals, na.rm = TRUE)
    z <- seq(min.raw, max.raw, length = length(col))
    #z <- quantile(vals, probs = seq(0,1,length = length(col)))
    par(mai=c(0.25,0.01,0.25,0.9))
    #Optimal margins depend upon graphics window
    #Try to exit if the plot cannot be produced, so that the device does not remain open
    #image(z = matrix(z, ncol = 1), col = col, breaks = seq(min.raw, max.raw, length.out=nlevels+1),
    #      xaxt = "n", yaxt = "n")
    #image(z = matrix(z, ncol = 1), col = col, breaks = quantile(vals, probs = seq(0, 1, length.out=nlevels+1)),
    #      xaxt = "n", yaxt = "n")
    #par(usr = c(0, 1, 0, 1))
    lv <- pretty(seq(min.raw, max.raw, length.out=nlevels+1))
    lv1 <- seq(min.raw, max.raw, length.out=nlevels+1)
    #print(lv)
    lv2 <- quantile(vals, probs = seq(0,1, length.out = nlevels+1), na.rm = T)
    #print(lv2)
    scale01 <- function(x, low = min(x), high = max(x)) {
      x <- (x - low)/(high - low)
      x
    }
    xv <- scale01(as.numeric(lv), min.raw, max.raw)
    xv1 <- scale01(as.numeric(lv1), min.raw, max.raw)
    image(z = matrix(z, nrow = 1), col = col, breaks = seq(min.raw, max.raw, length.out=nlevels+1),
          xaxt = "n", yaxt = "n")
    xvadj = 1/(2*nlevels)
    #axis(1, at = round(xv1,3), labels = round(lv2,3), cex.axis=0.8)
    #axis(4, at = seq(0-xvadj,1+xvadj,length.out=(nlevels+1)), labels = round(lv2,3), cex.axis=0.8,
    #     las=1)
    axis(4, at = seq(0-xvadj,1+xvadj,length.out=(nlevels+1)), labels = round(forbreaks,3),
         cex.axis=0.8, las=1)
    #axis(1, at = xv/2, labels = lv, cex.axis=0.8)
    mtext(colbarlab,at=0.5)
  }
  return(plotobj)
}

#' Plot a cladogram with RERs shown as either labels (type="label") or a color
#' heatmap along the branches (type="color")
#' Wraps around \code{\link{returnRersAsTree}} or \code{\link{treePlotNew}}, respectively
#'
#' @param treesObj. A treesObj created by \code{\link{readTrees}}
#' @param rermat. A residual matrix, output of the getAllResiduals() function
#' @param index. A character denoting the name of gene, or a numeric value
#' corresponding to the gene's row index in the residuals matrix
#' @param type. Whether to display RERs as branch labels ('label') or a heatmap ('color')
#' @param phenv. A phenotype vector returned by \code{\link{tree2Paths}} or \code{\link{foreground2Paths}}
#' @param .... Additional parameters to be passed to \code{\link{returnRersAsTree}}
#' or \code{\link{treePlotNew}}
#' #' @param figwid. Adjust x limits of plot.phylo by 1/figwid. May be related to figure width
#' (requires some optimization)
#' @return Plots a cladogram of the master tree with RERs displayed as branch labels or colors
#' @export

treePlotRers <- function(treesObj, rermat=NULL, index=NULL, type=c("label","color"),
                         phenv=NULL, figwid=10, ...) {
  type = match.arg(type)
  if (type=="label") {
    #use returnRersAsTree with plot = TRUE
    tmpout = returnRersAsTree(treesObj, rermat, index, phenv, plot = T, ...)
  }
  if (type=="color") {
    #use treePlotNew with edges colored by RER
    rerstoplot = returnRersAsTree(treesObj, rermat, index, phenv, plot = F)
    tmpout = treePlotNew(treesObj$trees[[index]], rerstoplot$edge.length,
                         colbarlab="RER", figwid=figwid, ...)
  }
  #return(tmpout)
}

#Plotting function using ggtree, with branch colors (currently) from edge lengths
#Demonstrates how to map edge order in R object to edge order in ggtree
treePlotGG = function(traittree, tiplabels = FALSE, title=NULL) {
  #use ggtree to make a pretty circular tree with branch lengths colored
  #currently bases colors on edge lengths in trait tree
  require(ggtree)
  forcol = c(rep("black",length(traittree$edge.length)))
  forcol[which(traittree$edge.length == -1)] = "turquoise"
  forcol[which(traittree$edge.length == 1)] = "red"
  forcol[which(is.na(traittree$edge.length))] = "gray"
  forcol = c(forcol,"goldenrod") #check to make sure the extra color is not included

  #Colors appear to be assigned in the following order:
  #1) terminal branches (in order of tip labels)
  #2) internal branches (by number of child node)
  #This means there may be a number skipped, and that they need to be mapped to edges somehow.
  #Make a vector of which color index corresponds to which edge index
  forcolgg = c(rep("black",max(traittree$edge)))
  for (g in c(1:length(forcolgg))) {
    #get the color based on the original forcol
    wgg = which(traittree$edge[,2] == g)
    if (length(wgg) > 0) {
      forcolgg[g] = forcol[wgg]
    }
  }

  #Title plots using ggtitle (optionally)
  treeplot = ggtree(traittree,layout='circular',branch.length='none',color=forcolgg)
  if (!is.null(title)) {
    treeplot = treeplot + ggtitle(title)
  }
  if (tiplabels) {
    treeplot = treeplot + geom_tiplab(aes(angle=angle),size=1.5)
  }
  return(treeplot)
}

#' Produce a gene tree with branch lengths representing RERs and (optionally)
#' display these RERs as branch labels
#'
#' @param treesObj. A treesObj created by \code{\link{readTrees}}
#' @param rermat. A residual matrix, output of the getAllResiduals() function
#' @param index. A character denoting the name of gene, or a numeric value corresponding to the gene's row index in the residuals matrix
#' @param phenv. A phenotype vector returned by \code{\link{tree2Paths}} or \code{\link{foreground2Paths}}
#' @param rer.cex. Numeric expansion for RER labels
#' @param tip.cex. Numeric expansion for tip labels
#' @param nalab. Label given to any NA RERs
#' @param plot. Whether to produce a plot displaying the RERs on the gene tree
#' @return An object of class "phylo" with edge lengths representing RERs for the given gene
#' @return If plot = TRUE, also displays a plot of the gene tree with edges labeled with RERs
#' @export

returnRersAsTree <- function(treesObj, rermat, index, phenv = NULL, rer.cex = 0.7,
                             tip.cex = 0.7, nalab = 'NA', plot = T){
  trgene <- treesObj$trees[[index]]
  ee=edgeIndexRelativeMaster(trgene, treesObj$masterTree)
  ii= treesObj$matIndex[ee[, c(2,1)]]
  rertree=unname(rermat[index,ii]) #avoid storing names
  rertree[is.nan(rertree)]=NA #replace NaNs from C functions
  if (plot) {
    trgene$edge.length <- rep(2,nrow(trgene$edge))
    par(mar = c(1,1,1,0))
    edgcols <- rep('black', nrow(trgene$edge))
    edgwds <- rep(1, nrow(trgene$edge))
    if(!is.null(phenv)){
      edgcols <- rep('black', nrow(trgene$edge))
      edgwds <- rep(1, nrow(trgene$edge))
      edgcols[phenv[ii]==1] <- 'red'
      edgwds[phenv[ii]==1] <- 2
    }
    plot.phylo(trgene, font = 2, edge.color = edgcols, edge.width = edgwds, cex = tip.cex)
    rerlab <- round(rertree,3)
    rerlab[is.na(rerlab)] <- nalab
    edgelabels(rerlab, bg = NULL, adj = c(0.5,0.9), col = edgcols, frame = 'none',cex = rer.cex, font =2)
  }
  trgene$edge.length <- rertree
  return(trgene)
}

#' Produce a vector of newick strings representing gene trees where the edge lengths
#' correspond to RERs
#'
#' @param treesObj. A treesObj created by \code{\link{readTrees}}
#' @param rermat. A residual matrix, output of the getAllResiduals() function
#' @return A named character vector of newick strings, one per gene, representing RERs as edge lengths'
#' @export
#'
returnRersAsNewickStrings <- function(treesObj, rermat){
  rerNwkstrings <- sapply(names(treesObj$trees), function(index){
    trgene <- treesObj$trees[[index]]
    ee=edgeIndexRelativeMaster(trgene, treesObj$masterTree)
    ii= treesObj$matIndex[ee[, c(2,1)]]
    rertree=rermat[index,ii]
    rertree[is.nan(rertree)]=NA
    trgene$edge.length <- rertree
    write.tree(trgene)
  })
}

#' Produce a multiPhylo object of all gene trees with branch lengths representing RERs
#' for each gene
#' @param treesObj. A treesObj created by \code{\link{readTrees}}
#' @param rermat. A residual matrix, output of the getAllResiduals() function
#' @return An object of class "multiPhylo" of named gene trees with edge lengths
#' representing RERs for the given gene
#' @export

returnRersAsTreesAll <- function(treesObj, rermat){
  whichgenes = intersect(names(treesObj$trees),rownames(rermat)) #allows subsetting
  allrers = lapply(whichgenes,returnRersAsTree,treesObj=treesObj,
                   rermat=rermat,plot=F)
  names(allrers)=whichgenes
  class(allrers)<-"multiPhylo"
  return(allrers)
}

#' Plot the residuals reflecting the relative evolutionary rates (RERs) of a gene across
#' species present in the gene tree
#' @param rermat. A residual matrix, output of the getAllResiduals() function
#' @param index. A character denoting the name of gene, or a numeric value corresponding to the gene's row index in the residuals matrix
#' @param phenv. A phenotype vector returned by \code{\link{tree2Paths}} or \code{\link{foreground2Paths}}
#' @return A plot of the RERs with foreground species labelled in red, and the rest in blue
#' @export

plotRers <- function(rermat=NULL, index= NULL, phenv = NULL, rers= NULL, method = 'k', xlims = NULL, plot = 1, xextend = 0.2, sortrers = F){
     if(is.null(rers)){
      e1 = rermat[index,][!is.na(rermat[index,])]
      colids = !is.na(rermat[index,])
      e1plot <- e1
      #print(e1plot)
      if(exists('speciesNames')){
          names(e1plot) <- speciesNames[names(e1),]
      }
      if(is.numeric(index)){
         gen = rownames(rermat)[index]
      }else{
         gen = index
      }
     }else{
      e1plot = rers
      gen = 'rates'
     }
     names(e1plot)[is.na(names(e1plot))]=""
     if(!is.null(phenv)){
      phenvid = phenv[colids]
      fgdcor = getAllCor(rermat[index,,drop=F],phenv, method = method)
      plottitle = paste0(gen, ': rho = ',round(fgdcor$Rho,4),', p = ',round(fgdcor$P,4))
      fgd = setdiff(names(e1plot)[phenvid == 1],"")
      df <- data.frame(species = names(e1plot), rer = e1plot, stringsAsFactors=FALSE) %>%
          mutate(mole = as.factor(ifelse(phenvid > 0,2,1)))
     }else{
      plottitle = gen
      fgd = NULL
      df <- data.frame(species = names(e1plot), rer = e1plot, stringsAsFactors=FALSE) %>%
          mutate(mole = as.factor(ifelse(0,2,1)))
     }
     #print(plottitle)
     if(sortrers){
      df = filter(df, species!="") %>%
          arrange(desc(rer))
     }
     #print(df)
     #df <- data.frame(species = names(e1plot), rer = e1plot, stringsAsFactors=FALSE) %>%
     #     mutate(mole = as.factor(ifelse(names(e1plot) %in% fgd,2,1)))
     #ll=c(min(df$rer)*1.1, max(df$rer)+xextend)
     if(is.null(xlims)){
          ll=c(min(df$rer)*1.1, max(df$rer)+xextend)
     }else{
          ll=xlims
     }
     g  <- ggplot(df, aes(x = rer, y=factor(species, levels = unique(ifelse(rep(sortrers, nrow(df)), species[order(rer)], sort(unique(species)))) ), col=mole, label=species)) + scale_size_manual(values=c(3,3))+ geom_point(aes(size=mole))+
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
     if(plot){
      print(g)
     }
     else{
      g
     }
}


nvmaster <- function(treesObj, useSpecies = NULL, fgd = NULL, plot = 0){
     #treesObj = simtrees
     #useSpecies = NULL
     #fgd = paste0('species',fgdbranchnums1)
     #fgd = matrix(paste0('species',fgdcomb),4)
     control = NULL
     if (is.null(useSpecies)){
          useSpecies=treesObj$masterTree$tip.label
     }
     both=intersect(treesObj$master$tip.label, useSpecies)
     allreport=treesObj$report[,both]
     ss=rowSums(allreport)
     iiboth=which(ss==length(both))
     ee=edgeIndexRelativeMaster(treesObj$masterTree, treesObj$masterTree)
     ii= treesObj$matIndex[ee[, c(2,1)]]
     allbranch=treesObj$paths[iiboth,ii]
     nv=t(projection(t(allbranch), method="AVE", returnNV = T))
     nv=as.vector(nv)
     mastertree=treesObj$master
     mastertree$edge.length=nv
     nn=character(length(nv))
     iim=match(1:length(treesObj$masterTree$tip.label), treesObj$masterTree$edge[,2])
     nn[iim]=treesObj$masterTree$tip.label
     names(mastertree$edge.length) = nn
     nvplot2 = sort(mastertree$edge.length)
     nvplot=nvplot2[names(nvplot2)!=""]
     if(plot){
          barcols = rep('black',length(nv))
          avlcols <- c('red','green','blue','yellow','orange','purple')
          nfgd = ifelse(!is.null(dim(fgd)),dim(fgd)[1],1)
          if(nfgd > 1){
               for(ii in 1:nfgd){
                    barcols[names(nvplot) %in% fgd[ii,]] = avlcols[ii]
               }
          }else{
               barcols[names(nvplot) %in% fgd] = avlcols[1]
          }
          bpl <- barplot(nvplot, horiz = T, col = barcols, names.arg = F, cex.axis = 1,
                         cex.lab = 2, xlab = 'average rate', space = 2, width = 0.5,
                         xlim = c(0,1.07*max(nvplot)))
          text( nvplot+rep(0.00,length(nvplot)) , bpl, labels = names(nvplot), srt = 0, pos = 4, cex = 0.75)
     }
     return(mastertree)
}


#' Plot the provided tree, (optionally) rerooted, with specified branches highlighted
#'
#' @param tree. A tree object.
#' @param outgroup. A vector of species to use to root the tree. If not provided, the tree remains unrooted.
#' @param hlspecies. A vector of species whose terminal branches to highlight, or a vector of branch numbers within the tree.
#' @param hlcols. Colors to use in highlighting the branches. If not specified, will use default R colors.
#' @param main. Main text for plot.
#' @return A plot of the the (optionally rerooted) tree, with branches highlighted.
#' @export

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
