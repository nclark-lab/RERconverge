CanonicalFormRenumberInternalNodes=function(tree){
     #par(mfrow=c(1,2))
     #  plot(tree)
     oo=order(tree$tip.label)
     tree$tip.label=tree$tip.label[oo]
     ii=match(1:length(oo), tree$edge[,2])
     tree$edge[ii,2]=order(oo)
     #plot(tree)
     renameintnodestoCanonical(rotateConstr(tree, sort(tree$tip.label)))
}
renameintnodestoCanonical <- function(x){
     #x <- CanonicalForm(t2)
     #tipmax=6
     tipmax <- length(x$tip)
     nodesinedge <- as.vector(t(x$edge))
     intids <- which(nodesinedge>tipmax)
     intnodesinedge <- nodesinedge[intids]
     intnodemin <- tipmax+1
     unqintnodesinedge <- intnodemin:(intnodemin+length(unique(intnodesinedge))-1)
     names(unqintnodesinedge) <- as.character(unique(intnodesinedge))
     newintnodesinedge <- unqintnodesinedge[as.character(intnodesinedge)]
     newnodesinedge <- nodesinedge
     newnodesinedge[intids] <- newintnodesinedge
     x$edge <- matrix(newnodesinedge, nrow = nrow(x$edge),byrow = T)
     x
}
CanonicalForm=function(tree){
     #par(mfrow=c(1,2))
     #  plot(tree)
     oo=order(tree$tip.label)
     tree$tip.label=tree$tip.label[oo]
     ii=match(1:length(oo), tree$edge[,2])
     tree$edge[ii,2]=order(oo)
     #plot(tree)
     rotateConstr(tree, sort(tree$tip.label))
}