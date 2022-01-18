#'Generates a binary phenotype tree using the list of tip foreground animals, foreground common ancestors, and their phylogenetic relationships
#' @param fg_vec A vector containing the tip foreground species
#' @param sisters_list A list containing pairs of "sister species" in the foreground set (put NULL if empty)
#' @param trees treesObj from \code{\link{readTrees}}
#' @param plotTree A boolean indicator for plotting the output tree (default=FALSE)
#' @param useSpecies An array containing the tip labels in the output tree
#' @return fg.tree A binary phenotype tree corresponding to the input information
#' @export
foreground2TreeClades=function(fg_vec,sisters_list=NULL,trees,plotTree=T,useSpecies=NULL){
  res.list = getForegroundInfoClades(fg_vec,sisters_list=sisters_list,trees,plotTree=plotTree,useSpecies=useSpecies)
  fg.tree = res.list$tree
  fg.tree
}

#'Generates a binary phenotype tree and foreground clades information using the list of tip foreground animals, the presence of foreground common ancestors, and their phylogenetic relationships
#' @param fg_vec A vector containing the tip foreground species
#' @param sisters_list A list containing pairs of "sister species" in the foreground set (put NULL if empty)
#' @param trees treesObj from \code{\link{readTrees}}
#' @param plotTree A boolean indicator for plotting the output tree (default=FALSE)
#' @param useSpecies An array containing the tip labels in the output tree
#' @return output.list A list containing 1) "tree" = a binary phenotype tree corresponding to the input information, 2) "fg.sisters.table" = a table containing all sister species in the foreground set
#' @export
getForegroundInfoClades=function(fg_vec,sisters_list=NULL,trees,plotTree=T,useSpecies=NULL){
  if (length(useSpecies)==0){
    useSpecies = trees$masterTree$tip.label
  }
  fg_tree = foreground2Tree(fg_vec,trees,plotTree=F,clade="all",useSpecies=useSpecies)
  edge = fg_tree$edge
  edge.length=fg_tree$edge.length

  ind.fg.edge = which(edge.length == 1)
  nodeIds.fg.edge = edge[ind.fg.edge,]

  tip.sisters = vector("integer",length=0)
  for (i in 1:length(sisters_list)){
    sisters = sisters_list[[i]]
    nodeId.sisters = which(useSpecies %in% sisters)
    if (length(nodeId.sisters)>0){
      tip.sisters = c(tip.sisters,nodeId.sisters)
    }
  }
  # Find and correct the pairs
  fg.sisters.table = matrix(nrow=0,ncol=2)
  colnames(fg.sisters.table) = c("species1","species2")

  if (length(as.vector(nodeIds.fg.edge)) > 2){
    nodeId.ca = sort(nodeIds.fg.edge[,1])
    fg_ca = vector("integer",length=0) # node IDs of the common ancestor foregrounds
    for (nn in 1:(length(nodeId.ca)-1)){
      if (nodeId.ca[nn] == nodeId.ca[nn+1]){
        nodeId.desc = nodeIds.fg.edge[which(nodeIds.fg.edge[,1]==nodeId.ca[nn]),2]
        if (length(which(nodeId.desc %in% tip.sisters)) > 0){
          fg_ca = c(fg_ca,nodeId.ca[nn])
          fg.sisters.table = rbind(fg.sisters.table, nodeIds.fg.edge[which(nodeIds.fg.edge[,1]==nodeId.ca[nn]),2])
        } else {
          fg_tree$edge.length[which(edge[,2] == nodeId.ca[nn])] = 0
        }
      }
    }
    rownames(fg.sisters.table) = fg_ca
    if (plotTree==T){
      plot(fg_tree)
    }
  }
  output.list = list("fg.sisters.table"=fg.sisters.table,"tree"=fg_tree)
  output.list
}

#'Calculates permuted correlation and enrichment statistics for binary phenotype
#' @param numperms An integer number of permulations
#' @param fg_vec A vector containing the tip foreground species
#' @param sisters_list  A list containing pairs of "sister species" in the foreground set (put NULL if empty)
#' @param root_sp The species to root the tree on
#' @param RERmat An RER matrix calculated using \code{\link{getAllResiduals}}.
#' @param trees treesObj from \code{\link{readTrees}}
#' @param mastertree A rooted, fully dichotomous tree derived from the treesObj master tree from \code{\link{readTrees}}.  Must not contain species not in traitvec
#' @param permmode Mode of binary permulation ("cc" for Complete Cases (default), "ssm" for Species Subset Match)
#' @param method statistical method to use for correlations (set to "k" (default) for Kendall Tau test)
#' @param min.pos minimum number of foreground species (default 2)
#' @param trees_list A list containing the trees of all genes of interest (formatted like trees in treesObj from \code{\link{readTrees}})
#' @param calculateenrich A boolean variable indicating if null permulation p-values for enrichment statistics
#' @param annotlist Pathway annotations
#' @return A list object with enrichment statistics, correlation p-val, rho, and correlation effect size
#' @export
getPermsBinary=function(numperms, fg_vec, sisters_list, root_sp, RERmat, trees, mastertree, permmode="cc", method="k", min.pos=2, trees_list=NULL, calculateenrich=F, annotlist=NULL){
  pathvec = foreground2Paths(fg_vec, trees, clade="all",plotTree=F)
  col_labels = colnames(trees$paths)
  names(pathvec) = col_labels

  if (permmode=="cc"){
    print("Running CC permulation")

    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhen(trees$masterTree, numperms, trees, root_sp, fg_vec, sisters_list, pathvec, permmode="cc")
    permulated.fg = mapply(getForegroundsFromBinaryTree, permulated.binphens[[1]])
    permulated.fg.list = as.list(data.frame(permulated.fg))
    phenvec.table = mapply(foreground2Paths,permulated.fg.list,MoreArgs=list(treesObj=trees,clade="all"))
    phenvec.list = lapply(seq_len(ncol(phenvec.table)), function(i) phenvec.table[,i])

    print("Calculating correlations")
    corMatList = lapply(phenvec.list, correlateWithBinaryPhenotype, RERmat=RERmat)

    #make enrich list/matrices to fill
    permPvals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permPvals)=rownames(RERmat)
    permRhovals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permRhovals)=rownames(RERmat)
    permStatvals=data.frame(matrix(ncol=numperms, nrow=nrow(RERmat)))
    rownames(permStatvals)=rownames(RERmat)

    for (i in 1:length(corMatList)){
      permPvals[,i] = corMatList[[i]]$P
      permRhovals[,i] = corMatList[[i]]$Rho
      permStatvals[,i] = sign(corMatList[[i]]$Rho)*-log10(corMatList[[i]]$P)
    }

  } else if (permmode=="ssm"){
    print("Running SSM permulation")

    if (is.null(trees_list)){
      trees_list = trees$trees
    }

    RERmat = RERmat[match(names(trees_list), rownames(RERmat)),]

    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhenSSMBatched(trees_list,numperms,trees,root_sp,fg_vec,sisters_list,pathvec)

    # Get species membership of the trees
    df.list = lapply(trees_list,getSpeciesMembershipStats,masterTree=masterTree,foregrounds=fg_vec)
    df.converted = data.frame(matrix(unlist(df.list), nrow=length(df.list), byrow=T),stringsAsFactors=FALSE)
    attr = attributes(df.list[[1]])
    col_names = attr$names
    attr2 = attributes(df.list)
    row_names = attr2$names

    colnames(df.converted) = col_names
    rownames(df.converted) = row_names

    df.converted$num.fg = as.integer(df.converted$num.fg)
    df.converted$num.spec = as.integer(df.converted$num.spec)

    spec.members = df.converted$spec.members

    # Group gene trees based on the similarity of their species membership
    grouped.trees = groupTrees(spec.members)
    ind.unique.trees = grouped.trees$ind.unique.trees
    ind.unique.trees = unlist(ind.unique.trees)
    ind.tree.groups = grouped.trees$ind.tree.groups

    # For each unique tree, produce a permuted tree. We already have this function, but we need a list of trees to feed in.
    unique.trees = trees_list[ind.unique.trees]

    # precompute clade mapping for each unique tree
    unique.map.list = mapply(matchAllNodesClades,unique.trees,MoreArgs=list(treesObj=trees))

    # calculate paths for each permulation
    unique.permulated.binphens = permulated.binphens[ind.unique.trees]
    unique.permulated.paths = calculatePermulatedPaths_apply(unique.permulated.binphens,unique.map.list,trees)

    permulated.paths = vector("list", length = length(trees_list))
    for (j in 1:length(permulated.paths)){
      permulated.paths[[j]] = vector("list",length=numperms)
    }
    for (i in 1:length(unique.permulated.paths)){
      ind.unique.tree = ind.unique.trees[i]
      ind.tree.group = ind.tree.groups[[i]]
      unique.path = unique.permulated.paths[[i]]
      for (k in 1:length(ind.tree.group)){
        permulated.paths[[ind.tree.group[k]]] = unique.path
      }
    }
    attributes(permulated.paths)$names = row_names

    print("Calculating correlations")
    RERmat.list = lapply(seq_len(nrow(RERmat[])), function(i) RERmat[i,])
    corMatList = mapply(calculateCorPermuted,permulated.paths,RERmat.list)
    permPvals = extractCorResults(corMatList,numperms,mode="P")
    rownames(permPvals) = names(trees_list)
    permRhovals = extractCorResults(corMatList,numperms,mode="Rho")
    rownames(permRhovals) = names(trees_list)
    permStatvals = sign(permRhovals)*-log10(permPvals)
    rownames(permStatvals) = names(trees_list)

  } else {
    stop("Invalid binary permulation mode.")
  }

  if (calculateenrich){
    realFgtree = foreground2TreeClades(fg_vec, sisters_list, trees, plotTree=F)
    realpaths = tree2PathsClades(realFgtree, trees)
    realresults = getAllCor(RERmat, realpaths, method=method, min.pos=min.pos)
    realstat =sign(realresults$Rho)*-log10(realresults$P)
    names(realstat) = rownames(RERmat)
    realenrich = fastwilcoxGMTall(na.omit(realstat), annotlist, outputGeneVals=F)

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

    counter=1;
    while (counter <= numperms){
      stat = permStatvals[,counter]
      names(stat) = rownames(RERmat)
      enrich=fastwilcoxGMTall(na.omit(stat), annotlist, outputGeneVals=F)
      #sort and store enrichment results
      groups=length(enrich)
      c=1
      while(c<=groups){
        current=enrich[[c]]
        enrich[[c]]=current[order(rownames(current)),]
        enrich[[c]]=enrich[[c]][match(rownames(permenrichP[[c]]), rownames(enrich[[c]])),]
        permenrichP[[c]][,counter]=enrich[[c]]$pval
        permenrichStat[[c]][,counter]=enrich[[c]]$stat
        c=c+1
      }
      counter = counter+1
    }
  }

  if(calculateenrich){
    data=vector("list", 5)
    data[[1]]=permPvals
    data[[2]]=permRhovals
    data[[3]]=permStatvals
    data[[4]]=permenrichP
    data[[5]]=permenrichStat
    names(data)=c("corP", "corRho", "corStat", "enrichP", "enrichStat")
  } else {
    data=vector("list", 3)
    data[[1]]=permPvals
    data[[2]]=permRhovals
    data[[3]]=permStatvals
    names(data)=c("corP", "corRho", "corStat")
  }
  data
}

#' @keywords internal
getForegroundsFromBinaryTree=function(tree){
  nameEdgesPerms.tree = nameEdgesPerms(tree)
  edge.length = as.logical(tree$edge.length)
  foregrounds = nameEdgesPerms.tree[edge.length]
  ind.tip = which(foregrounds != "")
  foregrounds = foregrounds[ind.tip]
  return(foregrounds)
}

#' @keywords internal
nameEdgesPerms=function(tree){
  if (is.null(tree$tip.label)) {
    nn = NULL
  } else {
    nn=character(nrow(tree$edge))
    iim=match(1:length(tree$tip.label), tree$edge[,2])
    nn[iim]=tree$tip.label
  }
  nn
}

#'Produces one CC binary permulation for a gene
#' @param trees treesObj from \code{\link{readTrees}}
#' @param mastertree A rooted, fully dichotomous tree derived from the treesObj master tree from \code{\link{readTrees}}.
#' @param root_sp The species to root the tree on
#' @param fg_vec a vector containing the foreground species
#' @param sisters_list  A list containing pairs of "sister species" in the foreground set (put NULL if empty)
#' @param pathvec A path vector generated from the real set of foreground animals
#' @param plotTreeBool Boolean indicator for plotting the output tree (default=FALSE)
#' @return A CC binary permulated tree
#' @export
simBinPhenoCC=function(trees, mastertree, root_sp, fg_vec, sisters_list=NULL, pathvec, plotTreeBool=F){
  tip.labels = mastertree$tip.label
  res = getForegroundInfoClades(fg_vec,sisters_list,trees,plotTree=F,useSpecies=tip.labels)
  fg_tree = res$tree
  fg.table = res$fg.sisters.table

  fgnum = length(which(fg_tree$edge.length == 1))
  internal = nrow(fg.table)
  tips=fgnum-internal # the number of tips

  num.tip.sisters.real = length(which(as.vector(fg.table) <= length(tip.labels)))

  top = NA
  num.tip.sisters.fake = 10000

  while(num.tip.sisters.fake!= num.tip.sisters.real){
    blsum=0
    while(blsum!=fgnum){
      t=root.phylo(trees$masterTree, root_sp, resolve.root = T)
      rm=ratematrix(t, pathvec)
      sims=sim.char(t, rm, nsim = 1)
      nam=rownames(sims)
      s=as.data.frame(sims)
      simulatedvec=s[,1]
      names(simulatedvec)=nam
      top=names(sort(simulatedvec, decreasing = TRUE))[1:tips]
      t=foreground2Tree(top, trees, clade="all", plotTree = F)
      blsum=sum(t$edge.length)

      t.table = findPairs(t)
      num.tip.sisters.fake = length(which(as.vector(t.table) <= length(tip.labels)))
    }
  }
  if (plotTreeBool){
    plot(t)
  }
  return(t)
}

#'Produces one SSM binary permulation for a gene
#' @param tree Tree of the gene of interest
#' @param trees treesObj from \code{\link{readTrees}}
#' @param root_sp The species to root the tree on
#' @param fg_vec A vector containing the foreground species
#' @param sisters_list A list containing pairs of "sister species" in the foreground set (put NULL if empty)
#' @param pathvec A path vector generated from the real set of foreground animals
#' @param plotTreeBool Boolean indicator for plotting the output tree (default=FALSE)
#' @return A SSM binary permulated tree
#' @export
simBinPhenoSSM=function(tree, trees, root_sp, fg_vec, sisters_list=NULL, pathvec, plotTreeBool=F){
  tip.labels = tree$tip.label # the set of species that exist in the gene tree
  ind_fg = which(tip.labels %in% fg_vec) # indices of the observed foreground animals that exist in the gene tree

  if (length(ind_fg) == 0){
    t = tree
    t$edge = NULL
    t$edge.length = NULL
    t$Nnode = NULL
    t$tip.label = NULL
  } else {
    fg_k = tip.labels[ind_fg] # the list of the observed foreground animals that exist in the gene tree

    res = getForegroundInfoClades(fg_k,sisters_list,trees,plotTree=F,useSpecies=tip.labels)
    fg_tree = res$tree
    fg.table = res$fg.sisters.table

    fgnum = length(which(fg_tree$edge.length == 1))
    internal = nrow(fg.table)
    tips=fgnum-internal # the number of tips

    num.tip.sisters.real = length(which(as.vector(fg.table) <= length(tip.labels)))

    top = NA
    num.tip.sisters.fake = 10000
    while(num.tip.sisters.fake!= num.tip.sisters.real){
      blsum=0
      while(blsum!=fgnum){
        t=root.phylo(trees$masterTree, root_sp, resolve.root = T) # roots the tree on the defined root
        rm=ratematrix(t, pathvec) # calculates the evolutionary VCV matrix
        sims=sim.char(t, rm, nsim = 1) # simulates evolution of discrete or continuous characters on a phylogenetic tree, results in phenotype vector
        nam=rownames(sims) # name of the animals
        s=as.data.frame(sims)
        simulatedvec=s[,1]
        names(simulatedvec)=nam
        top.all=names(sort(simulatedvec, decreasing = TRUE))
        top.tree_k = top.all[top.all %in% tip.labels]
        top = top.tree_k[1:tips]
        t=foreground2Tree(top, trees, clade="all", plotTree = F, useSpecies=tip.labels)
        blsum=sum(t$edge.length)

        t.table = findPairs(t)
        num.tip.sisters.fake = length(which(as.vector(t.table) <= length(tip.labels)))
      }
    }
  }
  if (plotTreeBool){
    plot(t)
  }
  return(t)
}

#' @keywords internal
findPairs=function(binary.tree){
  tip.labels = binary.tree$tip.label
  edge = binary.tree$edge
  edge.length = binary.tree$edge.length
  ind.fg.edge = which(edge.length == 1)
  fg.edges = edge[ind.fg.edge,]

  # Find the pairs
  fg.pairs.table = matrix(nrow=0,ncol=2)
  colnames(fg.pairs.table) = c("species1","species2")

  if (length(as.vector(fg.edges)) > 2){
    nodeId.start = sort(fg.edges[,1])
    fg_ca = vector("integer",length=0) # node IDs of the common ancestor foregrounds
    for (nn in 1:(length(nodeId.start)-1)){
      if (nodeId.start[nn] == nodeId.start[nn+1]){
        fg_ca = c(fg_ca,nodeId.start[nn])
        fg.pairs.table = rbind(fg.pairs.table, fg.edges[which(fg.edges[,1]==nodeId.start[nn]),2])
      }
    }
    rownames(fg.pairs.table) = fg_ca
  }
  fg.pairs.table
}

#'Produces binary permulations for a gene
#' @param tree Tree of the gene of interest (if permmode="cc", set this as the masterTree in trees (i.e., the output from \code{\link{readTrees}}))
#' @param numperms An integer number of permulations
#' @param trees treesObj from \code{\link{readTrees}}
#' @param root_sp The species to root the tree on
#' @param fg_vec A vector containing the foreground species
#' @param sisters_list  A list containing pairs of "sister species" in the foreground set (put NULL if empty)
#' @param pathvec A path vector generated from the real set of foreground animals
#' @return output.list a list containing the set of binary permulated trees
#' @export
generatePermulatedBinPhen=function(tree, numperms, trees, root_sp, fg_vec, sisters_list, pathvec, permmode="cc"){
  if (permmode=="cc"){
    tree_rep = lapply(1:numperms,rep_tree,tree=trees)
    permulated.binphens = lapply(tree_rep, simBinPhenoCC,mastertree=trees$masterTree,root_sp=root_sp, fg_vec=fg_vec,sisters_list=sisters_list,pathvec=pathvec,plotTreeBool=F)
  } else if (permmode=="ssm"){
    tree_rep = lapply(1:numperms,rep_tree,tree=tree)
    permulated.binphens = lapply(tree_rep,simBinPhenoSSM,trees=trees,root_sp=root_sp,fg_vec=fg_vec,sisters_list=sisters_list,pathvec=pathvec)
  } else {
    stop("Invalid binary permulation mode.")
  }
  output.list <- list()
  output.list[[1]] <- permulated.binphens
  return(output.list)
}
#' @keywords internal
rep_tree = function(num_input,tree){
  return(tree)
}

#'Produces binary SSM permulations for a list of genes
#' @param trees_list A list containing the trees of all genes of interest (formatted like trees in treesObj from \code{\link{readTrees}})
#' @param numperms An integer number of permulations
#' @param trees treesObj from \code{\link{readTrees}}
#' @param root_sp The species to root the tree on
#' @param fg_vec A vector containing the foreground species
#' @param sisters_list  A list containing pairs of "sister species" in the foreground set (put NULL if empty)
#' @param pathvec A path vector generated from the real set of foreground animals
#' @return simPhenoList A list containing binary permulated trees for each gene
#' @export
generatePermulatedBinPhenSSMBatched=function(trees_list,numperms,trees,root_sp,fg_vec,sisters_list,pathvec){
  masterTree = trees$masterTree
  master.tips = masterTree$tip.label
  df.list = lapply(trees_list,getSpeciesMembershipStats,masterTree=masterTree,foregrounds=fg_vec)
  df.converted = data.frame(matrix(unlist(df.list), nrow=length(df.list), byrow=T),stringsAsFactors=FALSE)
  attr = attributes(df.list[[1]])
  col_names = attr$names
  attr2 = attributes(df.list)
  row_names = attr2$names

  colnames(df.converted) = col_names
  rownames(df.converted) = row_names

  df.converted$num.fg = as.integer(df.converted$num.fg)
  df.converted$num.spec = as.integer(df.converted$num.spec)

  spec.members = df.converted$spec.members

  # Group gene trees based on the similarity of their species membership
  grouped.trees = groupTrees(spec.members)
  ind.unique.trees = grouped.trees$ind.unique.trees
  ind.unique.trees = unlist(ind.unique.trees)
  ind.tree.groups = grouped.trees$ind.tree.groups

  # For each unique tree, produce a permuted tree. We already have this function, but we need a list of trees to feed in.
  unique.trees = trees_list[ind.unique.trees]

  # Generate simulated phenotypes
  unique.pheno.list = mapply(generatePermulatedBinPhen,unique.trees,MoreArgs = list(numperms=numperms,trees=trees,root_sp=root_sp,fg_vec=fg_vec,sisters_list=sisters_list,pathvec=pathvec,permmode="ssm"))
  # Allocate the simulated phenotypes for unique trees to their respective groups
  simPhenoList = vector("list", length = length(trees_list))
  for (j in 1:length(simPhenoList)){
    simPhenoList[[j]] = vector("list",length=numperms)
  }
  for (i in 1:length(unique.pheno.list)){
    ind.unique.tree = ind.unique.trees[i]
    ind.tree.group = ind.tree.groups[[i]]
    unique.pheno = unique.pheno.list[[i]]
    for (k in 1:length(ind.tree.group)){
      simPhenoList[[ind.tree.group[k]]] = unique.pheno
    }
  }

  attributes(simPhenoList)$names = row_names

  return(simPhenoList)
}

#' @keywords internal
getSpeciesMembershipStats = function(tree,masterTree,foregrounds){
  master.tips = masterTree$tip.label
  tips = tree$tip.label
  spec_membership = which(master.tips %in% tips)
  fg_membership = which(foregrounds %in% tips)

  num_spec = length(spec_membership)
  num_fg = length(fg_membership)

  spec.members = rep(0,length(master.tips))
  spec.members[spec_membership] = 1
  spec.members = toString(spec.members)

  fg.members = rep(0,length(foregrounds))
  fg.members[fg_membership] = 1
  fg.members = toString(fg.members)

  df = data.frame("num.fg"=as.character(num_fg), "num.spec"=as.character(num_spec), "spec.members"=spec.members, "fg.members"=fg.members)

  return(df)
}

#' @keywords internal
groupTrees = function(spec.members){
  unique.trees = unique(spec.members)
  ind.tree.groups = lapply(unique.trees,findGroupedTrees,spec.members=spec.members)
  ind.unique.trees = lapply(ind.tree.groups, function(i) min(i))
  output.list = list()
  output.list$ind.unique.trees = ind.unique.trees
  output.list$ind.tree.groups = ind.tree.groups
  return(output.list)
}

#' @keywords internal
findGroupedTrees = function(unique.tree,spec.members){
  ind.grouped.trees = which(spec.members == unique.tree)
  return(ind.grouped.trees)
}

#'Calculates the clade mappings between the gene tree and the master tree (with the complete topology)
#' @param gene_tree A binary phenotype tree of a gene
#' @param treesObj treesObj from \code{\link{readTrees}}
#' @return output.map A list containing a dataframe of clades mapping
#' @export
matchAllNodesClades=function(gene_tree, treesObj){
  foregrounds = getForegroundsFromBinaryTree(gene_tree)
  tree1 = foreground2Tree(foregrounds,treesObj,clade="all",plotTree = F)

  map=matchNodesInject_mod(tree1,treesObj$masterTree)
  map=map[order(map[,1]),]
  #map

  output.map = list()
  output.map[[1]]=map
  output.map
}

#' @keywords internal
compareClades=function(clade2index,desc.tree2,clade1){
  output=NA
  clade2 = desc.tree2[[clade2index]]
  if (all(clade1 %in% clade2)){
    output = as.numeric(names(desc.tree2)[clade2index])
  }
  return(output)
}

#' @keywords internal
findMappedClade2Node=function(desc.tr1.index,desc.tr2.index.list,desc.tree1,desc.tree2){
  clade1 = desc.tree1[[desc.tr1.index]]
  mapped.clades.list = lapply(desc.tr2.index.list,compareClades,desc.tree2=desc.tree2,clade1=clade1)
  mapped.clades = unlist(mapped.clades.list)
  mapped.clades.nonNA = mapped.clades[!is.na(mapped.clades)]
  return(max(mapped.clades.nonNA))
}

#' @keywords  internal
matchNodesInject_mod=function (tr1, tr2){
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
  Nodes[,1] = as.numeric(names(desc.tr1))
  desc.tr1.index.list = as.list(1:length(desc.tr1))
  desc.tr2.index.list = as.list(1:length(desc.tr2))


  mapped.clade.list = lapply(desc.tr1.index.list,findMappedClade2Node,desc.tr2.index.list,desc.tree1=desc.tr1,desc.tree2=desc.tr2)
  Nodes[,2] = unlist(mapped.clade.list)

  iim=match(tr1$tip.label, tr2$tip.label)
  Nodes=rbind(cbind(1:length(tr1$tip.label),iim),Nodes)
  if(any(table(Nodes[,2])>1)){
    stop("Incorrect pseudorooting detected - use fixPseudoroot() function to correct trait tree topology")
  }

  Nodes
}

#'Calculates the paths for all permulated trees for a list of genes
#' @param permulated.trees.list A nested list of permulated phenotype trees for all the genes. The top layer lists the genes, and the nested layer lists the permulated binary trees (phylo objects) for each gene
#' @param map.list A list of maps corresponding to the genes listed in permulated.trees.list, the output of matchAllNodesClades
#' @param treesObj treesObj from \code{\link{readTrees}}
#' @return permulated.paths.list A nested list of permulated paths corresponding to the nested list of permulated phenotype trees
#' @export
calculatePermulatedPaths_apply=function(permulated.trees.list,map.list,treesObj){
  permulated.paths.list = mapply(calculatePermulatedPaths,permulated.trees.list,map.list,MoreArgs=list(treesObj=treesObj))
  permulated.paths.list
}

#' @keywords internal
calculatePermulatedPaths=function(permulated.trees,map,treesObj){
  permulated.paths=lapply(permulated.trees,tree2PathsClades,trees=treesObj)
  output.list = list()
  output.list[[1]] = permulated.paths
  output.list
}

#'A modification of the tree2Paths function that takes in pre-calculated mappings
#' @param tree the input tree to be converted into paths
#' @param trees treesObj from \code{\link{readTrees}}
#' @export
tree2PathsClades=function(tree,trees){
  map = matchAllNodesClades(tree,trees)
  path = tree2Paths_map(tree,map[[1]],trees)
  names(path) = colnames(trees$paths)
  path
}

#' @keywords internal
tree2Paths_map=function(tree, map, treesObj, binarize=NULL, useSpecies=NULL){
  if (class(tree)[1]=="phylo"){
    stopifnot(class(tree)[1]=="phylo")
    stopifnot(class(treesObj)[2]=="treesObj")

    if (is.null(tree$tip.label)){
      vals=as.double(rep(NA,length(treesObj$ap$dist)))
    } else {
      foregrounds = getForegroundsFromBinaryTree(tree)
      tree = foreground2Tree(foregrounds,treesObj,clade="all",plotTree = F)


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
      treePaths=allPaths(tree)

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
    }
  } else {
    vals=as.double(rep(NA,length(treesObj$ap$dist)))
  }

  vals
}

#' @keywords internal
allPaths=function(tree){
  dd=dist.nodes(tree) # returns a matrix with col_names and row_names denoting the numbers of the tips and the nodes, containing distances between nodes
  allD=double() # this is the 'path' vector that will be outputted in the end
  nn=matrix(nrow=0, ncol=2) # initialize an empty matrix that will store nodeIds of each node and its ancestors corresponding to each position in the path vector
  nA=length(tree$tip.label)+tree$Nnode # nA is the total number of nodes -- tree$Nnode is the number of internal nodes, length(tree$tip.label) is the number of tip nodes
  matIndex=matrix(nrow=nA, ncol=nA)
  index=1
  # Below here is where the vector of 'paths' are generated. Starting from the first to last tip labels, and then the internal nodes from root to last.
  for ( i in 1:nA){ # for each node i
    ia=getAncestors(tree,i) # find the ancestor nodes of node i
    if(length(ia)>0){ # If node i has ancestors at all
      allD=c(allD, dd[i, ia]) # extends allD per node, where each extension is the distances between node i and its ancestors
      nn=rbind(nn,cbind(rep(i, length(ia)), ia)) # extend the xxx-by-2 matrix by nodeId pairs of node i and its ancestors
      for (j in ia){
        matIndex[i,j]=index
        index=index+1
      }
    }
  }
  return(list(dist=allD, nodeId=nn, matIndex=matIndex))
}

#' @keywords internal
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

#'Calculate permulation correlation statistics
#' @param permulated.paths A nested list of permulated paths (e.g., output of \code{\link{calculatePermulatedPaths_apply}}
#' @param RERmat An RER matrix calculated using \code{\link{getAllResiduals}}.
#' @param min.sp Minimum number of species that must be present for a gene
#' @param min.pos Minimum number of species that must be present in the foreground (non-zero phenotype values)
#' @param method Method used to compute correlations. Accepts the same arguments as \code{\link{cor}}. Set to "auto" to select automatically based on the number of unique values in charP. This will also auto set the winsorization for Pearson correlation. Set winsorizetrait=some number and winsorizeRER=some number to override
#' @return A nested list containing the correlation statistics for the permulations
#' @export
calculateCorPermuted=function(permulated.paths,RERmat,min.sp=10,min.pos=2,method="k"){
  corMatList = lapply(permulated.paths,getAllCorSSM,RERmat,min.sp=min.sp,min.pos=min.pos,method=method)
  output.list <- list()
  output.list[[1]] <- corMatList
  return(output.list)
}

#' @keywords internal
getAllCorSSM=function(charP, RERmat, method="auto",min.sp=10, min.pos=2, winsorizeRER=NULL, winsorizetrait=NULL, weighted=F){
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
  dim(RERmat) <- c(1,length(RERmat))
  corout=matrix(nrow=nrow(RERmat), ncol=3)
  rownames(corout)=rownames(RERmat)

  colnames(corout)=c("Rho", "N", "P")

  for( i in 1:nrow(corout)){

    if(((nb<-sum(ii<-(!is.na(charP)&!is.na(RERmat[i,]))))>=min.sp)){
      if (method!="p"&&sum(charP[ii]!=0)<min.pos){
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

        cres=cor.test(x, y, method=method, exact=F)
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

#' @keywords internal
extractCorResults=function(corResultsList,numperms,mode="Rho"){
  output = matrix(NA,nrow=length(corResultsList),ncol=length(corResultsList[[1]]))
  for (i in 1:length(corResultsList)){
    gene = corResultsList[[i]]
    output[i,] = extractPermulationResults(gene,numperms,mode)
  }
  return(output)
}

#' @keywords internal
extractPermulationResults=function(gene,numperms,mode="Rho"){
  table_perm = lapply(gene,linearizeCorResults)
  df_perm = do.call(rbind,table_perm)
  if (mode=="Rho"){
    output=df_perm[,1]
  } else if (mode=="P"){
    output=df_perm[,3]
  }
  return(output)
}

#' @keywords internal
linearizeCorResults=function(cor_result){
  vec.cor = unlist(cor_result)
  return(vec.cor)
}

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
#' @param winR Integer winzorization value for RERmat
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
  names(realstat)=rownames(realresults)

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


  counter=1
  while(counter<=numperms){

    print(paste("running permutation: ", counter))

    #get correlation results
    out=getNullCor(traitvec, RERmat, mastertree, trees, type = type, winR=winR, winT=winT)
    stat=sign(out$Rho)*-log10(out$P)
    names(stat)=rownames(out)

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
        enrich[[c]]=enrich[[c]][match(rownames(permenrichP[[c]]), rownames(enrich[[c]])),]
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


#'Adaptively calculates permulation p-value for a single element
#' @param rer a row matrix for the given element (e.g., a row from the matrix output from \code{\link{getAllResiduals}})
#' @param permulated_foregrounds a list containing n sets of permulated foreground species names
#' @param observed_stats computed statistics for the observed trait (e.g., output from \code{\link{getAllCor}}, \code{\link{correlateWithBinaryPhenotype}}, or \code{\link{correlateWithContinuousPhenotype}})
#' @param alpha the significance level to control (default = 0.05)
#' @return permPval permulation (empirical) p-value
#' @export
adaptivePermulation=function(rer, permulated_foregrounds, observed_stats, alpha=0.05){
  dim(rer) = c(1, length(rer))
  observed_score = observed_stats$Rho

  max_permulations = length(permulated_foregrounds)
  maxnum_extreme = round(alpha*max_permulations)

  permulated_scores = rep(NA, length(permulated_foregrounds))

  for (k in 1:length(permulated_foregrounds)){
    #print(k)
    perm_path = foreground2Paths(permulated_foregrounds[k][[1]], trees, clade="terminal")
    perm_cor = correlateWithBinaryPhenotype(rer, perm_path)
    permulated_scores[k] = perm_cor$Rho
    computed_permulated_scores = permulated_scores[!is.na(permulated_scores)]

    if (length(computed_permulated_scores) >= 2* maxnum_extreme){
      median_null_scores = median(computed_permulated_scores)

      if (observed_score <= median_null_scores){
        one_sided_null_scores = computed_permulated_scores[which(computed_permulated_scores <= median_null_scores)]
        ind_extreme = which(one_sided_null_scores <= observed_score)
      } else if (observed_score > median_null_scores){
        one_sided_null_scores = computed_permulated_scores[which(computed_permulated_scores > median_null_scores)]
        ind_extreme = which(one_sided_null_scores >= observed_score)
      }

      if (length(ind_extreme) > maxnum_extreme || k == length(permulated_foregrounds)){
        permPval = min(maxnum_extreme+1, length(ind_extreme)+1)/(length(one_sided_null_scores)+1)
        break
      }
    }
  }
  permPval
}





#'Performs enrichment statistic permulations using existing gene correlation permulations
#' @param corperms Gene correlation permulations from \code{\link{getPermsContinuous}}
#' @param realenrich Pathway enrichment results using observed phenotype obtained from `correlateWithContinuousPhenotype` or `correlateWithBinaryPhenotype`
#' @param annotlist Pathway annotations
#' @return Full null permulation statistics and p-values for gene correlations and pathway enrichment
#' @export
getEnrichPerms=function(corperms, realenrich, annotlist){
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
      enrich[[c]]=enrich[[c]][match(rownames(permenrichP[[c]]), rownames(enrich[[c]])),]
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

#'Generates a permuted continuous phenotype given an observed continuous phenotype
#' @param namedvec A named numeric vector with phenotype values for each speices
#' @return A vector with permuted phenotype values
#' @export
permutevec=function(namedvec){
  #returns permuted vec
  n=names(namedvec)
  vec=sample(namedvec)
  names(vec)=n
  vec
}

#'Generates a simulated continuous phenotype given an observed continuous phenotype and a phylogeny
#' @param namedvec A named numeric vector with phenotype values for each speices
#' @param treewithbranchlengths A rooted phylogenetic tree with the same species as namedvec and branch lengths representing average evolutionary rate. The master tree from \code{\link{readTrees}} may be rooted and used for this parameter.
#' @return A vector with simulated phenotype values
#' @export
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

#'Generates a permulated continuous phenotype given an observed continuous phenotype and a phylogeny
#' @param namedvec A named numeric vector with phenotype values for each speices
#' @param treewithbranchlengths A rooted phylogenetic tree with the same species as namedvec and branch lengths representing average evolutionary rate. The master tree from \code{\link{readTrees}} may be rooted and used for this parameter.
#' @return A vector with permulated phenotype values
#' @export
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



#PGLS functions
####################

#'Generates a permulated phylogenetic tree with specified number of foreground branches.  User may specify the number of foreground branches that are internal branches.
#' @param trees treesObj output from \code{\link{readTrees}}
#' @param root Species on which to root the master tree
#' @param phenvec Named vector of 1's and 0's representing phenotype values for each species
#' @param fgnum Total number of foreground species - only required if internal foreground branches are required
#' @param internal Number of foreground species that should be internal branches - only required if internal foreground branches are required
#' @param drop Character vector (or single character variable) of species names to be removed from the master tree (such as species in trees but not in phenotype vector)
#' @return A tree with permulated phenotype values
#' @export
simBinPheno=function(trees, root, phenvec, fgnum=NULL, internal=0, drop=NULL){
  blsum=0
  if(is.null(fgnum)){
    fgnum=sum(phenvec)
  }
  tips=fgnum-internal
  while(blsum!=fgnum){
    t=root.phylo(trees$masterTree, root, resolve.root = T)
    t=drop.tip(t, drop)
    rm=ratematrix(t, phenvec)
    sims=sim.char(t, rm, nsim = 1)
    nam=rownames(sims)
    s=as.data.frame(sims)
    simulatedvec=s[,1]
    names(simulatedvec)=nam
    top=names(sort(simulatedvec, decreasing = TRUE))[1:tips]
    t=foreground2Tree(top, trees, clade="all", plotTree = F)
    blsum=sum(t$edge.length)
  }
  # plot(t)
  return(t)
}

#'Generates a permulated phenotype vector whose phylogeny matches a desired structure.  User may specify the number of foreground branches that are internal branches.
#' @param trees treesObj output from \code{\link{readTrees}}
#' @param root Species on which to root the master tree
#' @param phenvec Named vector of 1's and 0's representing phenotype values for each species
#' @param fgnum Total number of foreground species - only required if internal foreground branches are required
#' @param internal Number of foreground species that should be internal branches - only required if internal foreground branches are required
#' @param drop Character vector (or single character variable) of species names to be removed from the master tree (such as species in trees but not in phenotype vector)
#' @return A vector of permulated foreground species
#' @export
simBinPhenoVec=function(trees, root, phenvec, fgnum=NULL, internal=0, drop=NULL){
  blsum=0
  if(is.null(fgnum)){
    fgnum=sum(phenvec)
  }
  tips=fgnum-internal
  while(blsum!=fgnum){
    t=root.phylo(trees$masterTree, root, resolve.root = T)
    t=drop.tip(t, drop)
    rm=ratematrix(t, phenvec)
    sims=sim.char(t, rm, nsim = 1)
    nam=rownames(sims)
    s=as.data.frame(sims)
    simulatedvec=s[,1]
    names(simulatedvec)=nam
    top=names(sort(simulatedvec, decreasing = TRUE))[1:tips]
    t=foreground2Tree(top, trees, clade="all", plotTree = F)
    blsum=sum(t$edge.length)
  }
  # plot(t)
  return(top)
}



####################

#'Plots changes in number of statistically significant rate acceleration or deceleration versus the number of permulations
#' @param res correlation statistic output from \code{\link{correlateWithBinaryPhenotype}} or \code{\link{correlateWithContinuousPhenotype}}
#' @param perm.out output from \code{\link{getPermsBinary}} or \code{\link{getPermsContinuous}}
#' @param interval interval of number of permulations (e.g., interval = 10 means that number of positives with be calculated for number of permulations = 10, 20, 30, ...)
#' @param pvalthres p-value threshold for identifying statistically significant rate acceleration or deceleration
#' @param output.res Boolean defining whether to output a list object containing changes in the number of identified elements with changing number of permulations (default=FALSE)
#' @return A list containing changes in the number of identified elements with changing number of permulations (default=NULL)
#' @export
plotPositivesFromPermulations=function(res, perm.out, interval, pvalthres, output.res=FALSE){
  numperms_list = seq(interval, length(perm.out$corRho), interval)

  num_signif = NULL
  for (i in 1:length(numperms_list)){
    perm.i = list("corP"=perm.out$corP[,1:numperms_list[i]], "corRho"=perm.out$corRho[,1:numperms_list[i]], "corStat"=perm.out$corStat[,1:numperms_list[i]])
    permpval = permpvalcor(res, perm.i)
    num_signif = c(num_signif, length(which(permpval <= pvalthres)))
  }
  plot(numperms_list, num_signif, pch=19, xlab="no. of permulations", ylab="no. of significant elements", ylim=c(0, max(num_signif+1)))

  if (output.res == TRUE){
    out = list("num.permulations"=numperms_list, "num.significant"=num_signif)
  } else {
    out = NULL
  }
  out
}



