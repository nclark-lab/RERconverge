#'Produces the phenotype inputs for binary permulation functions
#' @param fgTree A binary phenotype tree containing the phenotype values (branch length = 1 for foreground species, 0 for background species)
#' @return out A list containing fg_vec and sisters_list for binary permulation functions
#' @export
getBinaryPermulationInputsFromTree=function(fgTree){
  unq_edge_lengths = unique(fgTree$edge.length)
  if (length(which(!(unq_edge_lengths %in% c(0,1)))) > 0){
    stop('Phenotype must be binary.')
  }
  all_edges = fgTree$edge
  num_tip_species = length(fgTree$tip.label)

  idx_fg_branches = which(fgTree$edge.length == 1)
  if (length(idx_fg_branches) == 1){
	fg_edges = fgTree$edge[idx_fg_branches,]
	fg_edges = t(as.data.frame(fg_edges))
	if (fg_edges[1,2] <=num_tip_species){
		tip_foregrounds = fgTree$tip.label[fg_edges[,2]]
	}
  } else {
	fg_edges = fgTree$edge[idx_fg_branches,]
	tip_fg_edges = fg_edges[which(fg_edges[,2] <= num_tip_species),]
	tip_foregrounds = fgTree$tip.label[tip_fg_edges[,2]]
  }
#  tip_fg_edges = fg_edges[which(fg_edges[,2] <= num_tip_species),]
#  tip_foregrounds = fgTree$tip.label[tip_fg_edges[,2]]

  idx_node_edges = which(fg_edges[,2] > num_tip_species)
  if (length(idx_node_edges) == 1){
    node_fg_edges = fg_edges[which(fg_edges[,2] > num_tip_species),]
    edge_i = node_fg_edges
    node_i = edge_i[2]
    # find daughters of node_i
    idx_daugthers_i = which(all_edges[,1] == node_i)
    daughter_nodeIds = all_edges[idx_daugthers_i,2]
    daughters = fgTree$tip.label[daughter_nodeIds]
    sisters_list = list('node_i'=daughters)
  } else if (length(idx_node_edges) == 0) {
    sisters_list = NULL
  } else {
    node_fg_edges = fg_edges[which(fg_edges[,2] > num_tip_species),]
    daughters_info_list = list()
    parents = NULL
    for (i in 1:nrow(node_fg_edges)){
      edge_i = node_fg_edges[i,]
      # find the daughters of this node
      idx_daughters_i = which(all_edges[,1] == edge_i[2])
      daughter_edges = all_edges[idx_daughters_i,]
      daughters_info_list[[i]] = daughter_edges[,2]
      parents = c(parents, edge_i[2])
    }
    names(daughters_info_list) = parents
    ### write something to order the branches based on depth
    tip_fg_ids = tip_fg_edges[,2]
    depth_order = rep(NA, length(daughters_info_list))
    names(depth_order) = names(daughters_info_list)
    order_assigned = NULL
    while(length(which(is.na(depth_order))) > 0){
      idx_na = which(is.na(depth_order))
      if (length(idx_na) > 0){
        for (j in 1:length(idx_na)){
          idx_na_j = idx_na[j]
          parent_j = parents[idx_na_j]
          daughters_j = daughters_info_list[[idx_na_j]]
          num_tip_daughters = length(which(daughters_j %in% tip_fg_ids))
          if (num_tip_daughters == 2){
            depth_order[idx_na_j] = 1
            order_assigned = c(order_assigned, parent_j)
          } else if (num_tip_daughters==1){
            node_daughter = daughters_j[which(daughters_j > length(fgTree$tip.label))]
            if (node_daughter %in% order_assigned){
              depth_order[idx_na_j] = 2
              order_assigned = c(order_assigned, parent_j)
            }
          } else if (num_tip_daughters==0){
            node_daughters = daughters_j[which(daughters_j > length(fgTree$tip.label))]
            if (length(which(node_daughters %in% order_assigned)) == 2){
              node_daughters_depths = depth_order[as.character(node_daughters)]
              depth_order[idx_na_j] = max(node_daughters_depths) + 1
              order_assigned = c(order_assigned, parent_j)
            }
          }
        }
      }
    }

    # construct the sisters list
    sisters_list = NULL
    counter=0
    unq_depth_order = sort(unique(depth_order))
    nodes_addressed = tip_fg_ids
    for (j in 1:length(unq_depth_order)){
      depth_order_j = depth_order[which(depth_order==unq_depth_order[j])]
      daughters_info_order_j = daughters_info_list[names(depth_order_j)]
      for (i in 1:length(daughters_info_order_j)){
        daughters_i = daughters_info_order_j[[i]]
        if (length(which(daughters_i <= length(fgTree$tip.label))) == 2){
          counter = counter+1
          tip_daughters = fgTree$tip.label[daughters_i]
          sisters_list[[counter]] = tip_daughters
          names(sisters_list)[counter] = names(daughters_info_order_j)[i]
          nodes_addressed = c(nodes_addressed, names(daughters_info_order_j)[i])
        } else if (length(which(daughters_i <=length(fgTree$tip.label))) == 1){
          counter = counter+1
          tip_daughter_id = daughters_i[which(daughters_i <= length(fgTree$tip.label))]
          tip_daughter = fgTree$tip.label[tip_daughter_id]
          node_daughter_id = daughters_i[which(daughters_i > length(fgTree$tip.label))]
          sisters_list[[counter]] = c(node_daughter_id, tip_daughter)
          names(sisters_list)[counter] = names(daughters_info_order_j)[i]
          nodes_addressed = c(nodes_addressed, names(daughters_info_order_j)[i])
        } else if (length(which(daughters_i <=length(fgTree$tip.label))) == 0){
          counter = counter+1
          sisters_list[[counter]] = daughters_i
          names(sisters_list)[counter] = names(daughters_info_order_j)[i]
          nodes_addressed = c(nodes_addressed, names(daughters_info_order_j)[i])
        }
      }
    }
  }
  out = list('fg_vec'=tip_foregrounds, 'sisters_list'=sisters_list)
  out
}


#'Generates a binary phenotype tree using the list of tip foreground animals, foreground common ancestors, and their phylogenetic relationships
#' @param fg_vec A vector containing the tip foreground species
#' @param sisters_list A list containing pairs of "sister species" in the foreground set (put NULL if empty)
#' @param trees treesObj from \code{\link{readTrees}}
#' @param plotTree A boolean indicator for plotting the output tree (default=FALSE)
#' @param useSpecies An array containing the tip labels in the output tree
#' @param transition A character string indicating whether transitions between background and foreground branches
#' are "bidirectional" or "unidirectional" (no foreground to background transitions, the default)
#' @return fg.tree A binary phenotype tree corresponding to the input information
#' @export
foreground2TreeClades=function(fg_vec,sisters_list=NULL,trees,plotTree=T,
                               useSpecies=NULL,transition="unidirectional"){
  res.list = getForegroundInfoClades(fg_vec,sisters_list=sisters_list,trees,
                                     plotTree=plotTree,useSpecies=useSpecies,
                                     transition=transition)
  fg.tree = res.list$tree
  fg.tree
}


#'Generates a binary phenotype tree and foreground clades information using the list of tip foreground animals, the presence of foreground common ancestors, and their phylogenetic relationships
#' @param fg_vec A vector containing the tip foreground species
#' @param sisters_list A list containing pairs of "sister species" in the foreground set (put NULL if empty)
#' @param trees treesObj from \code{\link{readTrees}}
#' @param plotTree A boolean indicator for plotting the output tree (default=FALSE)
#' @param useSpecies An array containing the tip labels in the output tree
#' @param transition A character string indicating whether transitions between background and foreground branches
#' are "bidirectional" or "unidirectional" (no foreground to background transitions, the default)
#' @return output.list A list containing 1) "tree" = a binary phenotype tree corresponding to the input information, 2) "fg.sisters.table" = a table containing all sister species in the foreground set
#' @export
getForegroundInfoClades=function(fg_vec,sisters_list=NULL,trees,plotTree=T,
                                 useSpecies=NULL,transition="unidirectional"){
  if (length(useSpecies)==0){
    useSpecies = trees$masterTree$tip.label
  }

  if (is.null(sisters_list)){
    fg.sisters.table=NULL
    fg_tree = foreground2Tree(fg_vec, trees, plotTree=F, clade="terminal",
                              useSpecies=useSpecies,transition=transition)
  } else {
    # Start with a temp phenotype tree assuming that all internal nodes are foregrounds
    fg_tree = foreground2Tree(fg_vec,trees,plotTree=F,clade="all",
                              useSpecies=useSpecies,transition=transition)
    edge = fg_tree$edge
    edge.length=fg_tree$edge.length

    ind.fg.edge = which(edge.length == 1)
    nodeIds.fg.edge = edge[ind.fg.edge,] # all foreground edges in the temp tree

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
      all.nodeId.ca = sort(nodeIds.fg.edge[,1])
      count_all_nodeId_ca = table(all.nodeId.ca)
      unq.nodeId.ca = unique(all.nodeId.ca)
      fg_ca = vector("integer",length=0) # node IDs of the common ancestor foregrounds
      nodes_addressed = NULL
      while (length(unq.nodeId.ca) != length(nodes_addressed)){
        nodeId.ca = sort(all.nodeId.ca[which(!(all.nodeId.ca %in% nodes_addressed))])
        if (length(nodeId.ca) == 1){
          nodes_addressed = c(nodes_addressed, nodeId.ca)
        } else {
          for (nn in 1:(length(nodeId.ca)-1)){

            if (nodeId.ca[nn] == nodeId.ca[nn+1]){
              nodeId.desc = nodeIds.fg.edge[which(nodeIds.fg.edge[,1]==nodeId.ca[nn]),2]
              if (length(which(nodeId.desc %in% tip.sisters)) > 0){
                fg_ca = c(fg_ca,nodeId.ca[nn])
                fg.sisters.table = rbind(fg.sisters.table, nodeIds.fg.edge[which(nodeIds.fg.edge[,1]==nodeId.ca[nn]),2])
                nodes_addressed = c(nodes_addressed, nodeId.ca[nn])
              } else {
                if (length(which(fg_tree$tip.label[nodeId.desc] %in% fg_vec)) == 2){
                  fg_tree$edge.length[which(edge[,2]==nodeId.ca[nn])] = 0
                  nodes_addressed = c(nodes_addressed, nodeId.ca[nn])
                } else {
                  if (length(which(nodeId.desc %in% nodes_addressed)) == 2){
                    fg_ca = c(fg_ca,nodeId.ca[nn])
                    fg.sisters.table = rbind(fg.sisters.table, nodeIds.fg.edge[which(nodeIds.fg.edge[,1]==nodeId.ca[nn]),2])
                    nodes_addressed = c(nodes_addressed, nodeId.ca[nn])
                  }
                }
              }
            } else {
              nodeId.desc = nodeIds.fg.edge[which(nodeIds.fg.edge[,1]==nodeId.ca[nn]),2]
              if (length(nodeId.desc) == 2){
                if (nodeId.ca[nn] != nodeId.ca[nn-1]){
                  fg_tree$edge.length[which(edge[,2] == nodeId.ca[nn])] = 0
                  nodes_addressed = c(nodes_addressed, nodeId.ca[nn])
                  nodes_addressed = unique(nodes_addressed)
                }
              } else {
                nodes_addressed = c(nodes_addressed, nodeId.ca[nn])
              }
            }
          }
        }
      }
      rownames(fg.sisters.table) = fg_ca
    }
  }
  if (plotTree==T){
    plot(fg_tree)
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
#' @param transition A character string indicating whether transitions between background and foreground branches
#' are "bidirectional" or "unidirectional" (no foreground to background transitions, the default)
#' @return A list object with enrichment statistics, correlation p-val, rho, and correlation effect size
#' @export
getPermsBinary=function(numperms, fg_vec, sisters_list, root_sp, RERmat, trees,
                        mastertree, permmode="cc", method="k", min.pos=2,
                        trees_list=NULL, calculateenrich=F, annotlist=NULL,
                        transition="unidirectional"){
  pathvec = foreground2Paths(fg_vec, trees, clade="all",plotTree=F,
                             transition=transition)
  col_labels = colnames(trees$paths)
  names(pathvec) = col_labels

  if (permmode=="cc"){
    print("Running CC permulation")

    print("Generating permulated trees")
    permulated.binphens = generatePermulatedBinPhen(trees$masterTree, numperms, trees, root_sp, fg_vec, sisters_list, pathvec, permmode="cc")
    permulated.fg = mapply(getForegroundsFromBinaryTree, permulated.binphens[[1]])
    permulated.fg.list = as.list(data.frame(permulated.fg))
    phenvec.table = mapply(foreground2Paths,permulated.fg.list,
                           MoreArgs=list(treesObj=trees,clade="all",transition=transition))
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
    df.list = lapply(trees_list,getSpeciesMembershipStats,masterTree=mastertree,foregrounds=fg_vec)
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
    unique.map.list = mapply(matchAllNodesClades,unique.trees,
                             MoreArgs=list(treesObj=trees,transition=transition))

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

#'Calculates permuted correlation and enrichment statistics for binary phenotype (does not enforce matching structure of sister species, but maintains internal & tips foreground counts within a fudge of the original tree)
#'@param fgdspecs A vector containing the tip foreground species
#'@param RERs An RER matrix calculated using \code{\link{getAllResiduals}}.
#'@param trees treesObj from \code{\link{readTrees}}
#'@param useSpecies A vector of species to include
#'@param ntrees An integer number of permulations
#'@param root The species to root the tree on
#'@param fudge The number up to which the permulated tree can differ in total foreground species
#'@param cors The gene correlation results from \code{\link{correlateWithBinaryPhenotype}}
#'@param phenvec A named vector of the phenotype
#'@param transition A character string indicating whether transitions between background and foreground branches
#' are "bidirectional" or "unidirectional" (no foreground to background transitions, the default)
#'@return A list object with enrichment statistics, correlation p-val, rho, and correlation effect size
#'@export
getPermsBinaryFudged <- function(fgdspecs, RERs, trees, useSpecies, ntrees, root, fudge = 5, cors,
                           phenvec, transition = "unidirectional") {
  # get counts on original tree
  t = foreground2Tree(fgdspecs, trees, plotTree = TRUE, clade = "all",
                      useSpecies = useSpecies, transition = transition)
  fgnum = sum(t$edge.length)
  tips = length(fgdspecs)
  internal = fgnum - tips

  # print summary
  print(paste("fgnum:", fgnum))
  print(paste("tips:", tips))
  print(paste("internal:", internal))

  # drop species in the tree that we don't want to use
  # this is the tree passed to the simulation function
  drop = trees$masterTree$tip.label[!(trees$masterTree$tip.label %in% useSpecies)]
  t=drop.tip(trees$masterTree, drop)
  t=root.phylo(t, root, resolve.root = T)

  # get the ratematrix
  rm=ratematrix(t, phenvec)

  # make the data frames
  statdf=data.frame(matrix(data=NA, nrow=nrow(RERs),ncol=ntrees),row.names=rownames(RERs))
  pvaldf=data.frame(matrix(data=NA, nrow=nrow(RERs),ncol=ntrees),row.names=rownames(RERs))

  # generate the trees
  count=1
  while(count<=ntrees){

    #get phenotype:
    blsum=0

    while(blsum>(fgnum+fudge) | blsum<(fgnum-fudge)){
      ###########################################
      sims=sim.char(t, rm, nsim = 1)[,,1] #sim.char returns a weird array data structure, [,,1] is the named vector we want
      top=names(sort(sims, decreasing = TRUE))[1:tips]
      ###########################################
      tf=foreground2Tree(top, trees, clade="all", plotTree = F, transition = transition)
      blsum=sum(tf$edge.length)
    }

    #get path:
    p=tree2Paths(tf, trees, useSpecies = useSpecies)

    #run correlation:
    c=correlateWithBinaryPhenotype(RERs, p)

    ###########################################
    # this assumes rownames will always match
    statdf[,count]=c$Rho
    pvaldf[,count]=c$P
    ###########################################

    print(paste0("finished perm: ", count))
    count=count+1
  }

  # get perm p-val:
  corswithpermp=cors
  rows=nrow(corswithpermp)
  corswithpermp$permP=rep(NA,rows)
  # this assumes rownames will always match
  for(g in 1:rows){
    if(is.na(corswithpermp$Rho[g])){
      p=NA
    }
    else{
      p=sum(abs(statdf[g,])>abs(cors$Rho[g]),na.rm=TRUE)/sum(!is.na(statdf[g,]))
    }
    corswithpermp$permP[g]=p
  }
  corswithpermp$permP.adj=p.adjust(corswithpermp$permP, method="BH")

  # return results
  return(list(res=corswithpermp, stat=statdf, pval=pvaldf))
}

#'Calculates permuted correlation and enrichment statistics for binary phenotype when performing an extant only analysis
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
#' @param transition A character string indicating whether transitions between background and foreground branches
#' are "bidirectional" or "unidirectional" (no foreground to background transitions, the default)
#' @return A list object with enrichment statistics, correlation p-val, rho, and correlation effect size
#' @export
getPermsBinaryExtantOnly=function(numperms, fg_vec, sisters_list, root_sp,
                                  RERmat, trees, mastertree, permmode="cc",
                                  method="k", min.pos=2, trees_list=NULL,
                                  calculateenrich=F, annotlist=NULL,
                                  transition="unidirectional"){
  pathvec = foreground2Paths(fg_vec, trees, clade="all",plotTree=F,
                             transition=transition)
  col_labels = colnames(trees$paths)
  names(pathvec) = col_labels

  print("Generating permulated trees")
  # list of permulated binary trees:
  permulated.binphens = generatePermulatedBinPhen(trees$masterTree, numperms, trees, root_sp, fg_vec, sisters_list, pathvec, permmode="cc")
  # matrix of fgd specs (each col is a list of fgd specs)
  permulated.fg = mapply(getForegroundsFromBinaryTree, permulated.binphens[[1]])
  # turn the matrix into a list (each entry is a set of fgd specs)
  permulated.fg.list = as.list(data.frame(permulated.fg))

  #phenvec.table = mapply(foreground2Paths,permulated.fg.list,MoreArgs=list(treesObj=trees,clade="all"))
  #phenvec.list = lapply(seq_len(ncol(phenvec.table)), function(i) phenvec.table[,i])

  # make a list of named phenotype vectors
  allspecs = na.omit(unique(col_labels))
  vec=rep(0, length(allspecs))
  names(vec) = allspecs

  phenvec.list = list()

  for(i in 1:length(permulated.fg.list)){
    phenvec.list[[i]] = vec
    phenvec.list[[i]][permulated.fg.list[[i]]] = 1
  }

  print("Calculating correlations")
  #corMatList = lapply(phenvec.list, correlateWithBinaryPhenotype, RERmat=RERmat)
  corMatList = lapply(phenvec.list, function(phen){getAllCorExtantOnly(RERmat, phen, method = "k", min.pos = min.pos)})

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

  if (calculateenrich){
    print("Calculating enrichments")
    #realFgtree = foreground2TreeClades(fg_vec, sisters_list, trees, plotTree=F)
    #realpaths = tree2PathsClades(realFgtree, trees)
    #realresults = getAllCor(RERmat, realpaths, method=method, min.pos=min.pos)
    realphen = vec
    realphen[fg_vec] = 1
    realresults = getAllCorExtantOnly(RERmat, realphen, method=method, min.pos=min.pos)
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
#' @param transition A character string indicating whether transitions between background and foreground branches
#' are "bidirectional" or "unidirectional" (no foreground to background transitions, the default)
#' @return A CC binary permulated tree
#' @export
simBinPhenoCC=function(trees, mastertree, root_sp, fg_vec, sisters_list=NULL,
                       pathvec, plotTreeBool=F, transition="unidirectional"){
  tip.labels = mastertree$tip.label
  res = getForegroundInfoClades(fg_vec,sisters_list,trees,plotTree=F,
                                useSpecies=tip.labels,transition=transition)  #### This function has problems
  fg_tree = res$tree
  fg.table = res$fg.sisters.table

  t=root.phylo(trees$masterTree, root_sp, resolve.root = T)
  rm=ratematrix(t, pathvec)

  if (!is.null(sisters_list)){
    fg_tree_info = getBinaryPermulationInputsFromTree(fg_tree)
    num_tip_sisters_true = unlist(fg_tree_info$sisters_list)
    num_tip_sisters_true = num_tip_sisters_true[which(num_tip_sisters_true %in% tip.labels)]
    num_tip_sisters_true = length(num_tip_sisters_true)
    fg_tree_depth_order = getDepthOrder(fg_tree)
  } else {
    fg_tree_depth_order = NULL
  }

  fgnum = length(which(fg_tree$edge.length == 1))
  if (!is.null(sisters_list)){
    internal = nrow(fg.table)
  } else {
    internal = 0
  }
  tips=fgnum-internal # the number of tips

  testcondition=FALSE
  while(!testcondition){
    blsum=0
    while(blsum!=fgnum){
      sims=sim.char(t, rm, nsim = 1)
      nam=rownames(sims)
      s=as.data.frame(sims)
      simulatedvec=s[,1]
      names(simulatedvec)=nam
      top=names(sort(simulatedvec, decreasing = TRUE))[1:tips]
      t_iter=foreground2Tree(top, trees, clade="all", plotTree = F,
                             transition = transition)
      blsum=sum(t_iter$edge.length)
    }
    t_info = getBinaryPermulationInputsFromTree(t_iter)
    if (!is.null(sisters_list)){
      num_tip_sisters_fake = unlist(t_info$sisters_list)
      num_tip_sisters_fake = num_tip_sisters_fake[which(num_tip_sisters_fake %in% tip.labels)]
      num_tip_sisters_fake = length(num_tip_sisters_fake)
      t_depth_order = getDepthOrder(t_iter)
      testcondition = setequal(sort(t_depth_order), sort(fg_tree_depth_order)) &&
        (num_tip_sisters_fake == num_tip_sisters_true)
    } else {
      t_depth_order = getDepthOrder(t_iter)
      testcondition = setequal(sort(t_depth_order), sort(fg_tree_depth_order))
    }
  }
  if (plotTreeBool){
    plot(t_iter)
  }
  return(t_iter)
}

#'Produces one CC binary permulation for a gene using midpoint rooting
#' @param trees treesObj from \code{\link{readTrees}}
#' @param mastertree A rooted, fully dichotomous tree derived from the treesObj master tree from \code{\link{readTrees}}.
#' @param fg_vec a vector containing the foreground species
#' @param sisters_list  A list containing pairs of "sister species" in the foreground set (put NULL if empty)
#' @param pathvec A path vector generated from the real set of foreground animals
#' @param plotTreeBool Boolean indicator for plotting the output tree (default=FALSE)
#' @param transition A character string indicating whether transitions between background and foreground branches
#' are "bidirectional" or "unidirectional" (no foreground to background transitions, the default)
#' @return A CC binary permulated tree
#' @export
simBinPhenoCCmidpoint=function(trees, mastertree, fg_vec, sisters_list=NULL,
                               pathvec, plotTreeBool=F, transition="unidirectional"){
  tip.labels = mastertree$tip.label
  res = getForegroundInfoClades(fg_vec,sisters_list,trees,plotTree=F,
                                useSpecies=tip.labels,transition=transition)  #### This function has problems
  fg_tree = res$tree
  fg.table = res$fg.sisters.table


  t = midpoint.root(mastertree)
  rm = ratematrix(t, pathvec)

  if (!is.null(sisters_list)){
    fg_tree_info = getBinaryPermulationInputsFromTree(fg_tree)
    num_tip_sisters_true = unlist(fg_tree_info$sisters_list)
    num_tip_sisters_true = num_tip_sisters_true[which(num_tip_sisters_true %in% tip.labels)]
    num_tip_sisters_true = length(num_tip_sisters_true)
    fg_tree_depth_order = getDepthOrder(fg_tree)
  } else {
    fg_tree_depth_order = NULL
  }

  fgnum = length(which(fg_tree$edge.length == 1))
  if (!is.null(sisters_list)){
    internal = nrow(fg.table)
  } else {
    internal = 0
  }
  tips=fgnum-internal # the number of tips

  testcondition=FALSE
  while(!testcondition){
    blsum=0
    while(blsum!=fgnum){
      sims=sim.char(t, rm, nsim = 1)
      nam=rownames(sims)
      s=as.data.frame(sims)
      simulatedvec=s[,1]
      names(simulatedvec)=nam
      top=names(sort(simulatedvec, decreasing = TRUE))[1:tips]
      t_iter=foreground2Tree(top, trees, clade="all", plotTree = F,
                             transition = transition)
      blsum=sum(t_iter$edge.length)
    }
    t_info = getBinaryPermulationInputsFromTree(t_iter)
    if (!is.null(sisters_list)){
      num_tip_sisters_fake = unlist(t_info$sisters_list)
      num_tip_sisters_fake = num_tip_sisters_fake[which(num_tip_sisters_fake %in% tip.labels)]
      num_tip_sisters_fake = length(num_tip_sisters_fake)
      t_depth_order = getDepthOrder(t_iter)
      testcondition = setequal(sort(t_depth_order), sort(fg_tree_depth_order)) &&
        (num_tip_sisters_fake == num_tip_sisters_true)
    } else {
      t_depth_order = getDepthOrder(t_iter)
      testcondition = setequal(sort(t_depth_order), sort(fg_tree_depth_order))
    }
  }
  if (plotTreeBool){
    plot(t_iter)
  }
  return(t_iter)
}



#' @keywords internal
getDepthOrder=function(fgTree){
  unq_edge_lengths = unique(fgTree$edge.length)
  if (length(which(!(unq_edge_lengths %in% c(0,1)))) > 0){
    stop('Phenotype must be binary.')
  }
  all_edges = fgTree$edge
  num_tip_species = length(fgTree$tip.label)

  idx_fg_branches = which(fgTree$edge.length == 1)
  if (length(idx_fg_branches)==1){
	fg_edges = fgTree$edge[idx_fg_branches,]
	fg_edges = t(as.data.frame(fg_edges))
  } else {
	fg_edges = fgTree$edge[idx_fg_branches,]
	tip_fg_edges = fg_edges[which(fg_edges[,2] <= num_tip_species),]
	tip_foregrounds = fgTree$tip.label[tip_fg_edges[,2]]
	node_fg_edges = fg_edges[which(fg_edges[,2] > num_tip_species),]
  }

  idx_node_edges = which(fg_edges[,2] > num_tip_species)
  if (length(idx_node_edges) == 1){
    node_fg_edges = fg_edges[which(fg_edges[,2] > num_tip_species),]
    node_fg_edges = t(as.data.frame(node_fg_edges))
  }
  if (length(idx_node_edges) == 0) {
    sisters_list = NULL
    depth_order = NULL
  } else {
    #node_fg_edges = fg_edges[which(fg_edges[,2] > num_tip_species),]
    daughters_info_list = list()
    parents = NULL
    for (i in 1:nrow(node_fg_edges)){
      edge_i = node_fg_edges[i,]
      # find the daughters of this node
      idx_daughters_i = which(all_edges[,1] == edge_i[2])
      daughter_edges = all_edges[idx_daughters_i,]
      daughters_info_list[[i]] = daughter_edges[,2]
      parents = c(parents, edge_i[2])
    }
    names(daughters_info_list) = parents
    ### write something to order the branches based on depth
    tip_fg_ids = tip_fg_edges[,2]
    depth_order = rep(NA, length(daughters_info_list))
    names(depth_order) = names(daughters_info_list)
    order_assigned = NULL
    while(length(which(is.na(depth_order))) > 0){
      idx_na = which(is.na(depth_order))
      if (length(idx_na) > 0){
        for (j in 1:length(idx_na)){
          idx_na_j = idx_na[j]
          parent_j = parents[idx_na_j]
          daughters_j = daughters_info_list[[idx_na_j]]
          num_tip_daughters = length(which(daughters_j %in% tip_fg_ids))
          if (num_tip_daughters == 2){
            depth_order[idx_na_j] = 1
            order_assigned = c(order_assigned, parent_j)
          } else if (num_tip_daughters==1){
            node_daughter = daughters_j[which(daughters_j > length(fgTree$tip.label))]
            if (node_daughter %in% order_assigned){
              depth_order[idx_na_j] = depth_order[as.character(node_daughter)] + 1
              order_assigned = c(order_assigned, parent_j)
            }
          } else if (num_tip_daughters==0){
            node_daughters = daughters_j[which(daughters_j > length(fgTree$tip.label))]
            if (length(which(node_daughters %in% order_assigned)) == 2){
              node_daughters_depths = depth_order[as.character(node_daughters)]
              depth_order[idx_na_j] = max(node_daughters_depths) + 1
              order_assigned = c(order_assigned, parent_j)
            }
          }
        }
      }
    }
  }
  depth_order
}



#'Produces one SSM binary permulation for a gene
#' @param tree Tree of the gene of interest
#' @param trees treesObj from \code{\link{readTrees}}
#' @param root_sp The species to root the tree on
#' @param fg_vec A vector containing the foreground species
#' @param sisters_list A list containing pairs of "sister species" in the foreground set (put NULL if empty)
#' @param pathvec A path vector generated from the real set of foreground animals
#' @param plotTreeBool Boolean indicator for plotting the output tree (default=FALSE)
#' @param transition A character string indicating whether transitions between background and foreground branches
#' are "bidirectional" or "unidirectional" (no foreground to background transitions, the default)
#' @return A SSM binary permulated tree
#' @export
simBinPhenoSSM=function(tree, trees, root_sp, fg_vec, sisters_list=NULL,
                        pathvec, plotTreeBool=F, transition="unidirectional"){
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

    res = getForegroundInfoClades(fg_k,sisters_list,trees,plotTree=F,
                                  useSpecies=tip.labels,transition=transition)
    fg_tree = res$tree
    fg.table = res$fg.sisters.table

    t=root.phylo(tree, root_sp, resolve.root = T)
    rm=ratematrix(t, pathvec)

    if (!is.null(sisters_list)){
      fg_tree_info = getBinaryPermulationInputsFromTree(fg_tree)
      num_tip_sisters_true = unlist(fg_tree_info$sisters_list)
      num_tip_sisters_true = num_tip_sisters_true[which(num_tip_sisters_true %in% tip.labels)]
      num_tip_sisters_true = length(num_tip_sisters_true)
      fg_tree_depth_order = getDepthOrder(fg_tree)
    }

    fgnum = length(which(fg_tree$edge.length == 1))
    if (!is.null(sisters_list)){
      internal = nrow(fg.table)
    } else {
      internal = 0
    }
    tips=fgnum-internal # the number of tips

    testcondition=FALSE
    while(!testcondition){
      blsum=0
      while(blsum!=fgnum){
        sims=sim.char(t, rm, nsim = 1)
        nam=rownames(sims)
        s=as.data.frame(sims)
        simulatedvec=s[,1]
        names(simulatedvec)=nam
        top.all=names(sort(simulatedvec, decreasing = TRUE))
        top.tree_k = top.all[top.all %in% tip.labels]
        top = top.tree_k[1:tips]
        t=foreground2Tree(top, trees, clade="all", plotTree = F,
                          useSpecies=tip.labels, transition=transition)
        blsum=sum(t$edge.length)
      }
      t_info = getBinaryPermulationInputsFromTree(t)
      if (!is.null(sisters_list)){
        num_tip_sisters_fake = unlist(t_info$sisters_list)
        num_tip_sisters_fake = num_tip_sisters_fake[which(num_tip_sisters_fake %in% tip.labels)]
        num_tip_sisters_fake = length(num_tip_sisters_fake)
        t_depth_order = getDepthOrder(t)
        testcondition = setequal(sort(t_depth_order), sort(fg_tree_depth_order)) &&
          (num_tip_sisters_fake == num_tip_sisters_true)
      } else {
        t_depth_order = getDepthOrder(t)
        testcondition = setequal(sort(t_depth_order), sort(fg_tree_depth_order))
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
#' @param transition A character string indicating whether transitions between background and foreground branches
#' are "bidirectional" or "unidirectional" (no foreground to background transitions, the default)
#' @return output.list a list containing the set of binary permulated trees
#' @export
generatePermulatedBinPhen=function(tree, numperms, trees, root_sp, fg_vec,
                                   sisters_list, pathvec, permmode="cc",
                                   transition="unidirectional"){
  if (permmode=="cc"){
    tree_rep = lapply(1:numperms,rep_tree,tree=trees)
    permulated.binphens = lapply(tree_rep, simBinPhenoCC,
                                 mastertree=trees$masterTree,root_sp=root_sp,
                                 fg_vec=fg_vec,sisters_list=sisters_list,
                                 pathvec=pathvec,plotTreeBool=F,transition=transition)
  } else if (permmode=="ssm"){
    tree_rep = lapply(1:numperms,rep_tree,tree=tree)
    permulated.binphens = lapply(tree_rep,simBinPhenoSSM,trees=trees,
                                 root_sp=root_sp,fg_vec=fg_vec,
                                 sisters_list=sisters_list,pathvec=pathvec,
                                 transition=transition)
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
#' @param transition A character string indicating whether transitions between background and foreground branches
#' are "bidirectional" or "unidirectional" (no foreground to background transitions, the default)
#' @return output.map A list containing a dataframe of clades mapping
#' @export
matchAllNodesClades=function(gene_tree, treesObj, transition="unidirectional"){
  foregrounds = getForegroundsFromBinaryTree(gene_tree)
  tree1 = foreground2Tree(foregrounds,treesObj,clade="all",plotTree = F,
                          transition=transition)

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
#' @param transition A character string indicating whether transitions between background and foreground branches
#' are "bidirectional" or "unidirectional" (no foreground to background transitions, the default)
#' @export
tree2PathsClades=function(tree,trees,transition="unidirectional"){
  map = matchAllNodesClades(tree,trees,transition=transition)
  path = tree2Paths_map(tree,map[[1]],trees,transition=transition)
  names(path) = colnames(trees$paths)
  path
}

#' @keywords internal
tree2Paths_map=function(tree, map, treesObj, binarize=NULL, useSpecies=NULL,
                        transition="unidirectional"){
  if (class(tree)[1]=="phylo"){
    stopifnot(class(tree)[1]=="phylo")
    stopifnot(class(treesObj)[2]=="treesObj")

    if (is.null(tree$tip.label)){
      vals=as.double(rep(NA,length(treesObj$ap$dist)))
    } else {
      foregrounds = getForegroundsFromBinaryTree(tree)
      tree = foreground2Tree(foregrounds,treesObj,clade="all",plotTree = F,
                             transition=transition)


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


#'Produces one binary permulation based on ranking of simulated branch lengths
#' @param trees treesObj from \code{\link{readTrees}}
#' @param root_sp The species to root the tree on
#' @param fg_vec a vector containing the foreground species
#' @param sisters_list  A list containing pairs of "sister species" in the foreground set (put NULL if empty)
#' @param plotTreeBool Boolean indicator for plotting the output tree (default=FALSE)
#' @param transition A character string indicating whether transitions between background and foreground branches
#' are "bidirectional" or "unidirectional" (no foreground to background transitions, the default)
#' @return A binary permulated tree
#' @export
simBinPhenoRank=function(trees, root_sp, fg_vec, sisters_list=NULL,
                         plotTreeBool=F, transition="unidirectional"){
  mastertree = trees$masterTree
  tip.labels = mastertree$tip.label
  res = getForegroundInfoClades(fg_vec,sisters_list,trees,plotTree=F,
                                useSpecies=tip.labels,transition=transition)
  fg_tree = res$tree
  pathvec = tree2PathsClades(fg_tree, trees)

  t=root.phylo(trees$masterTree, root_sp, resolve.root = T)
  ratem=ratematrix(t, pathvec)


  x = rnorm(n=length(t$edge.length))
  sd = sqrt(as.vector(ratem)*t$edge.length)

  y = matrix(0,nrow(t$edge),ncol(t$edge))
  alpha=0.1
  n = length(t$tip)
  for(i in 1:length(x)){
    if(t$edge[i,1]==(n+1))
      y[i,1]<-alpha # if at the root
    else
      y[i,1]<-y[match(t$edge[i,1],t$edge[,2]),2]
    y[i,2]<-y[i,1]+x[i]
  }

  rm(x)

  x<-c(y[1,1],y[,2])
  names(x)<-c(n+1,t$edge[,2])

  x<-x[as.character(1:(n+t$Nnode))]

  simphentree = t
  simedge = x[simphentree$edge[,2]]
  simphentree$edge.length = unname(simedge) + abs(min(simedge)) + 0.001

  numfg= sum(fg_tree$edge.length)

  simedgesort = sort(simphentree$edge.length, decreasing=T)
  simphenthreshold = simedgesort[numfg]

  bmphentree = simphentree
  simbinedge = simphentree$edge.length
  simbinedge[which(simbinedge < simphenthreshold)] = 0
  simbinedge[which(simbinedge >= simphenthreshold)] = 1
  bmphentree$edge.length = simbinedge

  if (plotTreeBool){
    plot(bmphentree)
  }

  return(bmphentree)
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

#'Calculates permuted correlation and enrichment statistics for an extant only analysis
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
getPermsContinuousExtantOnly=function(numperms, traitvec, RERmat, annotlist,
                                      trees, mastertree, calculateenrich=T,
                                      type="simperm", winR=3, winT=3,
                                      method="p", min.pos=0){

  #get real enrich and cors
  #realpaths=RERconverge::char2Paths(traitvec, trees)
  realresults=getAllCorExtantOnly(RERmat, traitvec, method = method, min.pos = min.pos, winsorizeRER = winR, winsorizetrait=winT)
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
    out=getNullCorExtantOnly(traitvec, RERmat, mastertree, trees, type = type, winR=winR, winT=winT)
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

#' @keywords internal
getNullCorExtantOnly=function(traitvec, RERmat, trimmedtree, genetrees, type="simperm", winR=NULL, winT=NULL){
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

  #get cor results
  # paths=RERconverge::char2Paths(na.omit(vec), genetrees)
  out=getAllCorExtantOnly(RERmat, na.omit(vec), method="p", min.pos=0, winsorizeRER=winR, winsorizetrait = winT)
  out
}

#'Adaptively calculates permulation p-value for a single element
#' @param rer a row matrix for the given element (e.g., a row from the matrix output from \code{\link{getAllResiduals}})
#' @param permulated_trees a list containing n sets of permulated trees (e.g., output from \code{\link{generatePermulatedBinPhen}})
#' @param observed_stats computed statistics for the observed trait (e.g., output from \code{\link{getAllCor}}, \code{\link{correlateWithBinaryPhenotype}}, or \code{\link{correlateWithContinuousPhenotype}})
#' @param alpha the significance level to control (default = 0.05)
#' @param mode "binary" for binary permulations, "continuous" for continuous permulations (default = "binary")
#' @return out data frame containing permPval (permulation p-value) and the corrected score (negative signifies deceleration, positive signifies acceleration)
#' @export
adaptivePermulation=function(rer, permulated_trees, observed_stats, alpha=0.05, mode="binary"){
  dim(rer) = c(1, length(rer))
  observed_score = observed_stats$Rho

  max_permulations = length(permulated_trees)
  maxnum_extreme = round(alpha*max_permulations)

  permulated_scores = rep(NA, length(permulated_trees))

  for (k in 1:length(permulated_trees)){
    perm_path = tree2Paths(permulated_trees[k][[1]], trees)
    if (mode == "binary"){
      perm_cor = correlateWithBinaryPhenotype(rer, perm_path)
    } else if (mode == "continuous"){
      perm_cor = correlateWithContinuousPhenotype(rer, perm_path)
    }
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

      if (length(ind_extreme) > maxnum_extreme || k == length(permulated_trees)){
        permPval = min(maxnum_extreme+1, length(ind_extreme)+1)/(length(one_sided_null_scores)+1)
        score = -log10(permPval)*sign(observed_score - median_null_scores)
        break
      }
    }
  }
  out = data.frame("permPval"=permPval, "score"=score)
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

#'Calculates correlation permutation pvals from output of \code{\link{getPermsContinuous}} and \code{\link{getPermsBinary}}
#' @param realcor Real enrichment statistics from \code{\link{fastwilcoxGMTall}}
#' @param permvals output from \code{\link{getPermsContinuous}} or \code{\link{getPermsBinary}}
#' @return A data frame containing permulation p-values and permulation statistics (positive denotes acceleration, negative denotes deceleration)
#' @export
permpvalcor=function(realcor, permvals){

  permcor=permvals$corRho
  realstat=realcor$Rho
  names(realstat)=rownames(realcor)

  permcor=permcor[match(names(realstat), rownames(permcor)),]

  permpvals=vector(length=length(realstat))
  names(permpvals)=names(realstat)
  permstats=vector(length=length(realstat))
  names(permstats)=names(realstat)
  count=1
  while(count<=length(realstat)){
    if(is.na(realstat[count])){
	  permpvals[count]=NA
    }else{
	  permcor_i = permcor[count,]
	  permcor_i = permcor_i[!is.na(permcor_i)]
	  if (length(permcor_i)==0){
	    permpvals[count]=NA
	    permstats[count]=NA
	  } else {
	    median_permcor = median(permcor_i)
	    if (realstat[count] >= median_permcor){
		  num = length(which(permcor_i >= realstat[count]))
		  denom = length(which(permcor_i >= median_permcor))
	    } else {
		  num = length(which(permcor_i <= realstat[count]))
		  denom = length(which(permcor_i <= median_permcor))
	    }
	    #num=sum(abs(permcor[count,])>abs(realstat[count]), na.rm=T)
	    #denom=sum(!is.na(permcor[count,]))
	    permpvals[count]=(num+1)/(denom+1)
	    permstats[count] = -log10(permpvals[count])*sign(realstat[count]-median_permcor)
	  }
    }
    count=count+1
  }
  permstats[which(is.na(permpvals))] = NA
  out = data.frame('permpval'=permpvals, 'permstats'=permstats)
  out
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
#' @param transition A character string indicating whether transitions between background and foreground branches
#' are "bidirectional" or "unidirectional" (no foreground to background transitions, the default)
#' @return A tree with permulated phenotype values
#' @export
simBinPheno=function(trees, root, phenvec, fgnum=NULL, internal=0, drop=NULL,
                     transition="unidirectional"){
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
    t=foreground2Tree(top, trees, clade="all", plotTree = F, transition=transition)
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
#' @param transition A character string indicating whether transitions between background and foreground branches
#' are "bidirectional" or "unidirectional" (no foreground to background transitions, the default)
#' @return A vector of permulated foreground species
#' @export
simBinPhenoVec=function(trees, root, phenvec, fgnum=NULL, internal=0, drop=NULL,
                        transition="unidirectional"){
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
    t=foreground2Tree(top, trees, clade="all", plotTree = F, transition=transition)
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

################################################################################
# Code for Permulations for Categorical Traits

# generates a set of N null tips with the number of species in each category matching the actual phenotype data
#' @keywords internal
getNullTips <- function(tree, Q, N, intlabels, root_prob = "stationary", percent_relax) {

  # GET TRUE TIP COUNTS
  true_counts = table(intlabels$mapped_states)

  # MAKE MATRIX TO STORE THE SETS OF NULL TIPS AND SETS OF INTERNAL NODES
  tips = matrix(nrow = N, ncol = length(tree$tip.label), dimnames = list(NULL, tree$tip.label))
  nodes = matrix(nrow = N, ncol = tree$Nnode)

  cnt = 0
  while(cnt < N) {
    # SIMULATE STATES
    sim = simulate_mk_model(tree, Q, root_probabilities = root_prob)
    sim_counts = table(sim$tip_states)

    # CHECK THAT ALL STATES GET SIMULATED IN THE TIPS
    if(length(unique(sim$tip_states)) < length(true_counts)) {
      next
    }

    # IF THE TIP COUNTS MATCH THE TIP COUNTS IN THE REAL DATA, ADD TO THE LIST
    # sum(true_counts == sim_counts) == length(true_counts)
    if(sum(abs(sim_counts - true_counts) <= true_counts*percent_relax) == length(true_counts)) {
      cnt = cnt + 1

      print(cnt)

      tips[cnt,] = sim$tip_states
      nodes[cnt,] = sim$node_states
    }
  }
  return(list(tips = tips, nodes = nodes))
}

# shuffles the categories around the tree based on ancestral likelihoods
# serves as a starting point for the function improveTree
#' @keywords internal
shuffleInternalNodes <- function(shuffled_states, available_nodes, ancliks, Nnode, ntips) {
  internal_states = vector(mode = "numeric", length = Nnode)

  if(length(shuffled_states) != nrow(ancliks)){
    stop("number of shuffled states and number of nodes with ancestral likelihoods do not match")
  }

  for(state in shuffled_states){
    if(length(available_nodes) > 1) {
      liks = ancliks[,state]
      node = sample(available_nodes, size = 1, prob = liks)
      available_nodes = available_nodes[- which(available_nodes == node)]
      ancliks = ancliks[- which(rownames(ancliks) == as.character(node)),]
      internal_states[node - ntips] = state
    }
    else { # only one node left
      internal_states[available_nodes - ntips] = state
    }
  }
  return(internal_states)
}

# For a set of null tips, shuffles the correct number of each category around the internal nodes
#' @keywords internal
getNullTrees <- function(node_states, null_tips, tree, Q) {

  nullTrees = list()

  for(i in 1:nrow(null_tips)) {
    print(i)
    tips = null_tips[i,]

    ancliks = getAncLiks(tree, tipvals = tips, Q = Q)

    ntips = length(tree$tip.label)
    available_nodes = (ntips + 1):(tree$Nnode + ntips)
    rownames(ancliks) = available_nodes

    shuffled_states = sample(node_states)
    internal_states = shuffleInternalNodes(shuffled_states,
                                           available_nodes = available_nodes,
                                           ancliks = ancliks,
                                           Nnode = tree$Nnode, ntips = ntips)
    tr = list(tips = tips, nodes = internal_states)
    nullTrees = append(nullTrees, list(tr))
  }
  return(nullTrees)
}

# rearranges the shuffled internal nodes to improve the likelihoods of the permulated trees
#' @keywords internal
improveTree <- function(tree, Q, P, nodes, tips, T0, Nk, cycles, alpha) {

  # get ancliks and max_states
  ancliks = getAncLiks(tree, tips, Q)

  states = c(tips, nodes)
  curr_lik = 1
  for(i in 1:nrow(tree$edge)){
    a = states[tree$edge[i,1]]
    d = states[tree$edge[i,2]]
    curr_lik = curr_lik * P[[i]][a, d]
  }

  # calculate initial ratios
  nstates = nrow(Q)
  ratios = c() # list of ratios
  ratio_info = matrix(nrow = (nstates - 1) * tree$Nnode, ncol = 3, dimnames = list(NULL, c("node", "state", "other.state"))) # info for each ratio

  # ns aren't the node numbers in the tree - they are the index of the internal node in nodes, node number in tree is n + ntips
  for(n in 1:tree$Nnode) {
    # calculate ratios
    pie = ancliks[n,]
    rr = pie[-nodes[n]] / pie[nodes[n]] # other states / state
    ratios = c(ratios, rr)
    # fill in ratio_info
    # rows = c((n-1)*3 + 1, (n-1)*3 + 2, (n-1)*3 + 3)
    rows = ((n-1)*(nstates-1) + 1):((n-1)*(nstates-1) + (nstates-1))
    ratio_info[rows,"node"] = rep(n, nstates - 1)
    ratio_info[rows,"state"] = rep(nodes[n], nstates - 1)
    ratio_info[rows,"other.state"] = (1:nstates)[-nodes[n]]
  }

  # pre-calculate and store edge numbers for each node
  ntips = length(tree$tip.label)
  edg_nums = lapply(seq_along(vector(mode = "list", length = tree$Nnode + ntips)), function(x){
    c(which(tree$edge[,1] == x),(which(tree$edge[,2] == x)))
  })

  j = 1 # iteration counter
  k = 1 # cycle counter
  Tk = T0

  while(k <= cycles) {

    # get 2 nodes to swap
    nn = nodes

    # 1: pick a node randomly, weighted by the ratios
    r1 = sample(1:length(ratios), 1, prob = ratios)
    n1 = ratio_info[r1, "node"] # node 1
    s1 = ratio_info[r1, "state"] # state1
    s2 = ratio_info[r1, "other.state"] # state2

    # 2: pick a node to swap it with
    ii = intersect(which(ratio_info[,"state"] == s2), which(ratio_info[,"other.state"] == s1))
    if(length(ii) > 1) {
      n2 = sample(ratio_info[ii,"node"], 1, prob = ratios[ii]) # node2
    } else { # only one node with state2
      n2 = ratio_info[ii,"node"]
    }

    # make the swap
    nn[n1] = s2
    nn[n2] = s1

    # calculate new likelihood
    states_new = c(tips, nn)
    states_old = c(tips, nodes)

    r = 1
    for(i in unique(c(edg_nums[[n1 + ntips]], edg_nums[[n2 + ntips]]))){ # check this over many cases including when n1 and n2 effect the same edge
      ao = states_old[tree$edge[i,1]]
      do = states_old[tree$edge[i,2]]

      an = states_new[tree$edge[i,1]]
      dn = states_new[tree$edge[i,2]]
      r = r * (P[[i]][an,dn] / P[[i]][ao, do])
    }

    if(r >= 1) { # if the swap increases likelihood, commit to the swap

      nodes = nn

      curr_lik = curr_lik * r # this should do the same thing, BUT CHECK THIS GETS THE SAME RESULT IN MULTIPLE CASES!

      # update ratios
      # rows1 = c((n1-1)*3 + 1, (n1-1)*3 + 2, (n1-1)*3 + 3) # rows to update ratios for n1
      rows1 = ((n1-1)*(nstates-1) + 1):((n1-1)*(nstates-1) + (nstates-1))

      pie = ancliks[n1,]
      rr = pie[-s2] / pie[s2] # other states / state
      ratios[rows1] = rr

      # fill in ratio_info
      ratio_info[rows1,"state"] = rep(s2, nstates - 1)
      ratio_info[rows1,"other.state"] = (1:nstates)[-s2]

      # rows2 = c((n2-1)*3 + 1, (n2-1)*3 + 2, (n2-1)*3 + 3) # rows to update ratios for n2
      rows2 = ((n2-1)*(nstates-1) + 1):((n2-1)*(nstates-1) + (nstates-1))

      pie = ancliks[n2,]
      rr = pie[-s1] / pie[s1] # other states / state
      ratios[rows2] = rr

      # fill in ratio_info
      ratio_info[rows2,"state"] = rep(s1, nstates - 1)
      ratio_info[rows2,"other.state"] = (1:nstates)[-s1]

    }
    else { # make jump with probability u

      # calculate u which includes dividing by tmp
      dh = -log(curr_lik * r) + log(curr_lik)
      u = exp(-dh/Tk)
      if(u == 0) warning("u is zero")

      # if(u == 0) stop(paste("temp is", Tk))

      if(runif(1) <= u) {

        nodes = nn

        curr_lik = curr_lik * r # CHECK THIS GETS THE SAME RESULT

        # update ratios
        # rows1 = c((n1-1)*3 + 1, (n1-1)*3 + 2, (n1-1)*3 + 3) # rows to update ratios for n1
        rows1 = ((n1-1)*(nstates-1) + 1):((n1-1)*(nstates-1) + (nstates-1))

        pie = ancliks[n1,]
        rr = pie[-s2] / pie[s2] # other states / state
        ratios[rows1] = rr
        # fill in ratio_info

        ratio_info[rows1,"state"] = rep(s2, nstates - 1)
        ratio_info[rows1,"other.state"] = (1:nstates)[-s2]

        # rows2 = c((n2-1)*3 + 1, (n2-1)*3 + 2, (n2-1)*3 + 3) # rows to update ratios for n2
        rows2 = ((n2-1)*(nstates-1) + 1):((n2-1)*(nstates-1) + (nstates-1))

        pie = ancliks[n2,]
        rr = pie[-s1] / pie[s1] # other states / state
        ratios[rows2] = rr
        # fill in ratio_info
        ratio_info[rows2,"state"] = rep(s1, nstates - 1)
        ratio_info[rows2,"other.state"] = (1:nstates)[-s1]
      }
    }

    # increment j
    j = j + 1

    # print(curr_lik)

    # move to next cycle if necessary
    if(j >= Nk) {
      j = 1 # reset j
      # Tk = T0 * alpha^k
      # Tk = T0 / (1 + alpha * log(k))
      Tk = T0 / (1 + alpha*k)
      k = k + 1
    }

  }
  end = Sys.time()
  return(list(nodes = nodes, lik = log10(curr_lik)))
}

#' @param treesObj trees object returned by readTrees
#' @param phenvals the named phenotype vector
#' @param rm the rate model, it should be the same as the one used to reconstruct the ancestral history of the trait
#' @param rp root prior, the default is auto
#' @param ntrees the number of null trees to generate
#' @param root_prob The probabilities of the different states at the root for the simulations. Can be "flat", "stationary", or a numeric vector of length Nstates. See the root_probabilities parameter under simulate_mk_model in the castor package for more information.
#' @return a set of permulated phenotype trees
#' @export
categoricalPermulations <- function(treesObj, phenvals, rm, rp = "auto", ntrees, percent_relax = 0){

  # check percent_relax is one value or a vector of length = # traits
  if(!(length(percent_relax) == 1 || length(percent_relax) == length(unique(phenvals)))) {
    stop("percent_relax is the wrong length")
  }

  # PRUNE TREE, ORDER PHENVALS, MAP TO STATE SPACE
  tree = treesObj$masterTree
  keep = intersect(names(phenvals), tree$tip.label)
  tree = pruneTree(tree, keep)
  phenvals = phenvals[tree$tip.label]
  intlabels = map_to_state_space(phenvals)

  # FIT A TRANSITION MATRIX ON THE DATA
  message("Fitting transition matrix")
  Q = fit_mk(tree, intlabels$Nstates, intlabels$mapped_states,
             rate_model = rm, root_prior = rp)$transition_matrix

  # GET NULL TIPS (AND STORE INTERNAL NODES FROM SIMULATIONS TOO)
  message("Simulating trees")
  simulations = getNullTips(tree, Q, ntrees, intlabels,
                            percent_relax = percent_relax)

  ancliks = getAncLiks(tree, intlabels$mapped_states, Q = Q)
  node_states = getStatesAtNodes(ancliks)

  # GET SHUFFLED STARTING-POINT TREES
  message("Shuffling internal states")
  nullTrees = getNullTrees(node_states, simulations$tips, tree, Q)

  P = lapply(tree$edge.length, function(x){expm(Q * x)})

  # IMPROVE LIKELIHOOD OF EACH NULL TREE
  message("Improving tree likelihoods")
  improvedNullTrees = lapply(nullTrees, function(x){
    list(tips = x$tips, nodes = improveTree(tree, Q, P, x$nodes, x$tips, 10, 10, 100, 0.9)$nodes)
  })

  # RETURN
  message("Done")
  return(list(sims = simulations, trees = improvedNullTrees, startingTrees = nullTrees))
}

#' @param realCors the output of correlateWithCategoricalPhenotype
#' @param nullPhens the list item named trees in the output of categoricalPermulations. If categoricalPermulations is run with extantOnly = TRUE, nullPhens is a matrix in which each row is a set of tips and getPermPvalsCategorical should be run with extantOnly = TRUE as well.
#' @param phenvals the named phenotype vector
#' @param treesObj the trees object returned by readTrees
#' @param RERmat the matrix of RERs returned by getAllResiduals, should be the same one used to calculate realCors
#' @param method either "kw" for Kruskal Wallis, the default, or "aov" for ANOVA
#' @return Permulation p-values for a categorical phenotype
#' @export
getPermPvalsCategorical <- function(realCors, nullPhens, phenvals,
                                    treesObj, RERmat, method = "kw",
                                    min.sp = 10, min.pos = 2,
                                    winsorizeRER=NULL, winsorizetrait=NULL,
                                    weighted=F,
                                    extantOnly = FALSE) {
  # CHECK IF TRAIT IS BINARY
  binary = FALSE
  if(method != "kw" & method != "aov"){
    message("Binary method provided. Setting binary to TRUE. Note: binary phenotype values should be TRUE and FALSE for correct results.")
    binary = TRUE
  }

  # PRUNE TREE
  tree = treesObj$masterTree
  keep = intersect(names(phenvals), tree$tip.label)
  tree = pruneTree(tree, keep)

  # UNROOT THE TREE IF IT IS ROOTED
  if (is.rooted(tree)) {
    tree = unroot(tree)
  }

  # GENERATE PATHS IF NOT EXTANT ONLY
  if(!extantOnly){
    # generate the paths
    if(!binary) {
      message("Generating null paths")
      nullPaths = lapply(nullPhens, function(x){
        tr = tree # make a copy of the tree
        tr$edge.length = c(x$tips, x$nodes)[tr$edge[,2]] # assign states to edges
        tree2Paths(tr, treesObj, categorical = TRUE, useSpecies = names(phenvals)) # calculate paths
      })
    } else {
      message("Generating null paths")
      nullPaths = lapply(nullPhens, function(x){
        tr = tree # make a copy of the tree
        tr$edge.length = c(x$tips, x$nodes)[tr$edge[,2]] # assign states to edges
        # subtract 1 to convert 1 - FALSE, 2 - TRUE to 0 and 1
        # THIS IS A QUICK FIX, WON'T WORK IF THE FOREGROUND IS 1 AND BACKGROUND IS 2!!!
        tree2Paths(tr, treesObj, categorical = TRUE, useSpecies = names(phenvals)) - 1 # calculate paths
      })
    }
  }

  # calculate correlation statistics
  message("Calculating correlation statistics")

  # make matrices to store the results
  if(!extantOnly){
    corsMatPvals = matrix(nrow = nrow(RERmat), ncol = length(nullPhens),
                          dimnames = list(rownames(RERmat), NULL))
    corsMatEffSize = matrix(nrow = nrow(RERmat), ncol = length(nullPhens),
                            dimnames = list(rownames(RERmat), NULL))
  } else {
    corsMatPvals = matrix(nrow = nrow(RERmat), ncol = nrow(nullPhens),
                          dimnames = list(rownames(RERmat), NULL))
    corsMatEffSize = matrix(nrow = nrow(RERmat), ncol = nrow(nullPhens),
                            dimnames = list(rownames(RERmat), NULL))
  }


  if(!binary){
    # make matrices for the pairwise testing
    Ppvals = lapply(1:length(realCors[[2]]), matrix, data = NA, nrow = nrow(RERmat),
                    ncol = length(nullPhens), dimnames = list(rownames(RERmat), NULL))
    names(Ppvals) = names(realCors[[2]])

    Peffsize = lapply(1:length(realCors[[2]]), matrix, data = NA, nrow = nrow(RERmat),
                      ncol = length(nullPhens), dimnames = list(rownames(RERmat), NULL))
    names(Peffsize) = names(realCors[[2]])
  }

  # run getAllCor on every null phenotype and store the pval/effect size for every gene
  if(extantOnly) {
    if(!binary) {
      for(i in 1:nrow(nullPhens)) {
        # phenvals is row i of nullPhens
        cors = getAllCorExtantOnly(RERmat, nullPhens[i,], method = method,
                                   min.sp = min.sp, min.pos = min.pos, winsorizeRER = winsorizeRER,
                                   winsorizetrait = winsorizetrait)
        corsMatPvals[,i] = cors[[1]]$P # store p values
        corsMatEffSize[,i] = cors[[1]]$Rho # store effect size

        # add results of pairwise tests
        for(j in 1:length(cors[[2]])){ # loop through each table in cors[[2]]
          Ppvals[[names(cors[[2]])[j]]][,i] = cors[[2]][[j]]$P
          Peffsize[[names(cors[[2]])[j]]][,i] = cors[[2]][[j]]$Rho
        }
      }
    } else {
      for(i in 1:nrow(nullPhens)) {
        print(i)
        # phenvals is row i of nullPhens
        cors = getAllCorExtantOnly(RERmat, nullPhens[i,], method = method,
                                   min.sp = min.sp, min.pos = min.pos, winsorizeRER = winsorizeRER,
                                   winsorizetrait = winsorizetrait)
        corsMatPvals[,i] = cors$P # store p values
        corsMatEffSize[,i] = cors$Rho # store effect size
      }
    }
  } else {
    if(!binary) {
      for(i in 1:length(nullPaths)) {
        cors = getAllCor(RERmat, nullPaths[[i]], method = method, min.sp=min.sp,
                         min.pos=min.pos, winsorizeRER=winsorizeRER,
                         winsorizetrait=winsorizetrait, weighted=weighted)
        corsMatPvals[,i] = cors[[1]]$P # store p values
        corsMatEffSize[,i] = cors[[1]]$Rho # store effect size

        # add results of pairwise tests
        for(j in 1:length(cors[[2]])){ # loop through each table in cors[[2]]
          Ppvals[[names(cors[[2]])[j]]][,i] = cors[[2]][[j]]$P
          Peffsize[[names(cors[[2]])[j]]][,i] = cors[[2]][[j]]$Rho
        }
      }
    } else {
      for(i in 1:length(nullPaths)) {
        cors = getAllCor(RERmat, nullPaths[[i]], method = method, min.sp=min.sp,
                         min.pos=min.pos, winsorizeRER=winsorizeRER,
                         winsorizetrait=winsorizetrait, weighted=weighted)
        corsMatPvals[,i] = cors$P # store p values
        corsMatEffSize[,i] = cors$Rho # store effect size
      }
    }
  }

  message("Obtaining permulations p-values")
  if(!binary){
    # calculate empirical pvals
    N = nrow(realCors[[1]])

    realCors[[1]]$permP = rep(NA, N) # add a column to the real results for the empirical pvals

    for(j in 1:length(realCors[[2]])){ # loop through each table in realCors[[2]]
      realCors[[2]][[j]]$permP = rep(NA, N)
    }

    for(gene in 1:N) {
      # check whether the gene is NA in realCors and if so, set p = NA
      if(is.na(realCors[[1]]$Rho[gene])) {
        p = NA
      } else {
        # count number of times null is more extreme than observed
        # effect size for Kruskal-Wallis is epsilon squared, only assumes non-negative values
        # ANOVA is still using eta2 (as of right now)
        p = sum(corsMatEffSize[gene,] > realCors[[1]]$Rho[gene], na.rm = TRUE) / sum(!is.na(corsMatEffSize[gene,]))
      }
      realCors[[1]]$permP[gene] = p

      for(j in 1:length(realCors[[2]])) {
        if(is.na(realCors[[2]][[j]]$Rho[gene])) {
          p = NA
        } else {
          # I think the effect size for Tukey and Dunn are bidirectional - using absolute value, NEED TO CHECK WITH AMANDA
          p = sum(abs(Peffsize[[names(realCors[[2]][j])]][gene,]) > abs(realCors[[2]][[j]]$Rho[gene]), na.rm = TRUE) / sum(!is.na(Peffsize[[names(realCors[[2]][j])]][gene,]))
        }
        realCors[[2]][[j]]$permP[gene] = p
      }
    }
  } else {
    # calculate empirical pvals
    N = nrow(realCors)

    realCors$permP = rep(NA, N) # add a column to the real results for the empircal pvals

    for(gene in 1:N) {
      if(is.na(realCors$Rho[gene])){
        p = NA
      } else {
        # count number of times null is more extreme than observed
        # use abs because stat can be positive or negative
        p = sum(abs(corsMatEffSize[gene,]) > abs(realCors$Rho[gene]), na.rm = TRUE) / sum(!is.na(corsMatEffSize[gene,]))
      }
      realCors$permP[gene] = p
    }
  }

  message("Done")
  # return results
  if(!binary){
    return(list(res = realCors, pvals = list(corsMatPvals,Ppvals), effsize = list(corsMatEffSize,Peffsize)))
  }
  else {
    return(list(res = realCors, pvals = corsMatPvals, effsize = corsMatEffSize))
  }
}

#' @param cors the output of correlateWithCategoricalPhenotype
#' @param annotlist the list of annotations from which to calculate enrichment statistics
#' @param outputGeneVals a boolean indicating whether to include gene names in the output
#' @return enrichment statistics for the categorical test and each pairwise test
#' @export
getRealEnrichments <- function(cors, annotlist, outputGeneVals = FALSE){
  # make a list to store all the enrichments
  enrich_list = list()

  # calculate enrichments for cors[[1]]
  enrich_list[[1]] = fastwilcoxGMTall(getStat(cors[[1]]), annotList = annotlist, outputGeneVals = outputGeneVals, alternative = "greater")

  # loop through and calculate enrichments for cors[[2]]
  enrich_list[[2]] = vector(mode = "list", length = length(cors[[2]]))
  names(enrich_list[[2]]) = names(cors[[2]])

  for(i in 1:length(cors[[2]])) {
    enrich_list[[2]][[i]] = fastwilcoxGMTall(getStat(cors[[2]][[i]]), annotList = annotlist, outputGeneVals = outputGeneVals)
  }
  # return the list of enrichments (match same kind of format as cors for consistency)
  return(enrich_list)
}

#' @param perms the output of getPermPvalsCategorical
#' @param realenrich the output of getRealEnrichments
#' @param annotlist the list of annotations from which to calculate enrichment statistics
#' @return enrichment statistics
#' @export
getEnrichPermsCategorical <- function(perms, realenrich, annotlist){
  res = vector(mode = "list", length = 2)

  # get enrichment results for KW/ANOVA on all categories
  res[[1]] = getEnrichAllCategories(perms, realenrich, annotlist)

  res[[2]] = vector(mode = "list", length = length(realenrich[[2]]))
  names(res[[2]]) = names(realenrich[[2]])

  # get enrichment results for each pairwise test
  for(name in names(res[[2]])){
    res[[2]][[name]] = getEnrichFromPairwiseTest(perms, realenrich, annotlist, name)
  }
  return(res)
}

# helper function that gets enrichment statistics for the KW/ANOVA results on all the categories
#' @keywords internal
getEnrichAllCategories <- function(perms, realenrich, annotlist){
  numperms = ncol(perms$pvals[[1]])

  # make the lists and data frames to store the results

  enrichP = vector(mode = "list", length = length(realenrich[[1]]))
  enrichStat = vector(mode = "list", length = length(realenrich[[1]]))

  names(enrichP) = names(realenrich[[1]])
  names(enrichStat) = names(realenrich[[1]])

  # give the rows in these matrices the same names as in the realenrich
  for(c in 1:length(realenrich[[1]])){
    enrichP[[c]] = data.frame(matrix(ncol = numperms, nrow = nrow(realenrich[[1]][[c]])))
    rownames(enrichP[[c]]) = rownames(realenrich[[1]][[c]])
    enrichStat[[c]] = data.frame(matrix(ncol = numperms, nrow = nrow(realenrich[[1]][[c]])))
    rownames(enrichStat[[c]]) = rownames(realenrich[[1]][[c]])
  }

  # calculate enrichment stats for each permulation
  for(count in 1:numperms){
    # print(count)
    # get the p value and enrich stat vectors
    P = perms$pvals[[1]][,count]
    effsize = perms$effsize[[1]][,count]

    # make a stat vector
    stat = getStat(data.frame(P = P, Rho = effsize))

    # run the enrichment
    enrich = fastwilcoxGMTall(stat, annotlist, outputGeneVals = FALSE, alternative="greater")

    # store the enrichment values
    for(c in 1:length(enrich)){
      #                                   the row indices in enrich[[c]] of the row names in enrichP
      enrichP[[c]][,count] = enrich[[c]][match(rownames(enrichP[[c]]), rownames(enrich[[c]])),]$pval
      enrichStat[[c]][,count] = enrich[[c]][match(rownames(enrichP[[c]]), rownames(enrich[[c]])),]$stat
    }
  }

  # return the results
  return(list(enrichP = enrichP, enrichStat = enrichStat))
}

# this function returns results for ONE pairwise test, it gets called mulitple times in the main function
# name is the name of the pairwise test e.g. "1 - 3"
#' @keywords internal
getEnrichFromPairwiseTest <- function(perms, realenrich, annotlist, name){
  numperms = ncol(perms$pvals[[1]])
  # make the lists and data frames to store the results

  enrichP = vector(mode = "list", length = length(realenrich[[2]][[name]]))
  enrichStat = vector(mode = "list", length = length(realenrich[[2]][[name]]))

  names(enrichP) = names(realenrich[[2]][[name]])
  names(enrichStat) = names(realenrich[[2]][[name]])

  # give the rows in these matrices the same names as in the realenrich
  for(c in 1:length(realenrich[[1]])){
    enrichP[[c]] = data.frame(matrix(ncol = numperms, nrow = nrow(realenrich[[2]][[name]][[c]])))
    rownames(enrichP[[c]]) = rownames(realenrich[[2]][[name]][[c]])
    enrichStat[[c]] = data.frame(matrix(ncol = numperms, nrow = nrow(realenrich[[2]][[name]][[c]])))
    rownames(enrichStat[[c]]) = rownames(realenrich[[2]][[name]][[c]])
  }

  # calculate enrichment stats for each permulation
  for(count in 1:numperms){
    print(count)
    # get the p value and enrich stat vectors
    P = perms$pvals[[2]][[name]][,count]
    effsize = perms$effsize[[2]][[name]][,count]

    # make a stat vector
    stat = getStat(data.frame(P = P, Rho = effsize))

    # run the enrichment
    enrich = fastwilcoxGMTall(stat, annotlist, outputGeneVals = FALSE)

    # store the enrichment values
    for(c in 1:length(enrich)){
      #                                   the row indices in enrich[[c]] of the row names in enrichP
      enrichP[[c]][,count] = enrich[[c]][match(rownames(enrichP[[c]]), rownames(enrich[[c]])),]$pval
      enrichStat[[c]][,count] = enrich[[c]][match(rownames(enrichP[[c]]), rownames(enrich[[c]])),]$stat
    }
  }

  # return the results
  return(list(enrichP = enrichP, enrichStat = enrichStat))
}

#' @param permenrich output of getEnrichPermsCategorical
#' @param realenrich output of getRealEnrichments
#' @param binary whether the trait is binary i.e. has two categories or not
#' @return the p-values for each enrichment category
#' @export
getEnrichPermPvals <- function(permenrich, realenrich, binary = FALSE){
  if(!binary){ # code for categorical traits with pairwise tests
    # calculate pvals for the enrichment of all categories

    # make a list of groups to store the pvalues
    pval_groups = vector(mode = "list", length = length(realenrich[[1]]))
    names(pval_groups) = names(realenrich[[1]])

    # loop through the groups of the enrichment
    for(c in 1:length(realenrich[[1]])){
      # make a list to store the pvals for that group
      pvals = c()
      for(i in 1:nrow(realenrich[[1]][[c]])){
        # for debugging - check that row names match, this should never happen because they were given the same row names
        if(rownames(realenrich[[1]][[c]])[i] != rownames(permenrich[[1]]$enrichStat[[c]])[i]){
          warning("row names between real enrichment and perm enrich stats do not match!")
        }
        # if the stat is NA, add NA to pval list
        if(is.na(realenrich[[1]][[c]]$stat[i])) {
          pvals = c(pvals, NA)
        }
        # otherwise calculate the p-value
        # do NOT use absolute value because it is a one-sided test
        else {
          p = sum(permenrich[[1]]$enrichStat[[c]][i,] > realenrich[[1]][[c]]$stat[i], na.rm = TRUE)
          p = p/sum(!is.na(permenrich[[1]]$enrichStat[[c]][i,]))
          pvals = c(pvals, p)
        }
      }
      # store the pvals in pval_groups list
      names(pvals) = rownames(realenrich[[1]][[c]])
      pval_groups[[c]] = pvals
    }

    # calculate pvals for the pairwise tests
    pw_pvals = vector(mode = "list", length = length(realenrich[[2]]))
    names(pw_pvals) = names(realenrich[[2]])

    for(n in 1:length(pw_pvals)){
      # make a list of groups to store the pvalues
      pval_groupstmp = vector(mode = "list", length = length(realenrich[[2]][[n]]))
      names(pval_groupstmp) = names(realenrich[[2]][[n]])

      # loop through the groups of the enrichment
      for(c in 1:length(realenrich[[2]][[n]])){
        # make a list to store the pvals for that group
        pvals = c()
        for(i in 1:nrow(realenrich[[2]][[n]][[c]])){
          # for debugging - check that row names match, this should never happen because they were given the same row names
          if(rownames(realenrich[[2]][[n]][[c]])[i] != rownames(permenrich[[2]][[names(realenrich[[2]])[n]]]$enrichStat[[c]])[i]){
            warning("row names between real enrichment and perm enrich stats do not match!")
          }
          # if the stat is NA, add NA to pval list
          if(is.na(realenrich[[2]][[n]][[c]]$stat[i])) {
            pvals = c(pvals, NA)
          }
          # otherwise calculate the p-value
          # use absolute value because it is a two-sided test
          else {
            p = sum(abs(permenrich[[2]][[names(realenrich[[2]])[n]]]$enrichStat[[c]][i,]) > abs(realenrich[[2]][[n]][[c]]$stat[i]), na.rm = TRUE)
            p = p/sum(!is.na(permenrich[[2]][[names(realenrich[[2]])[n]]]$enrichStat[[c]][i,]))
            pvals = c(pvals, p)
          }
        }
        # store the pvals in pval_groups list
        names(pvals) = rownames(realenrich[[2]][[n]][[c]])
        pval_groupstmp[[c]] = pvals
      }
      # store pval_groups in the pw_pvals list
      pw_pvals[[n]] = pval_groupstmp
    }

    return(list(pval_groups, pw_pvals))

  } else { # code for binary traits (permenrich is a list of lenght 2 with enrichP and enrichStat)
    pval_groups = vector(mode = "list", length = length(realenrich))
    names(pval_groups) = names(realenrich)

    # loop through the groups of the enrichment
    for(c in 1:length(realenrich)){
      # make a list to store the pvals for that group
      pvals = c()
      for(i in 1:nrow(realenrich[[c]])){
        # for debugging - check that row names match, this should never happen because they were given the same row names
        if(rownames(realenrich[[c]])[i] != rownames(permenrich$enrichStat[[c]])[i]){
          warning("row names between real enrichment and perm enrich stats do not match!")
        }
        # if the stat is NA, add NA to pval list
        if(is.na(realenrich[[c]]$stat[i])) {
          pvals = c(pvals, NA)
        }
        # otherwise calculate the p-value
        else {
          p = sum(abs(permenrich$enrichStat[[c]][i,]) > abs(realenrich[[c]]$stat[i]), na.rm = TRUE)
          p = p/sum(!is.na(permenrich$enrichStat[[c]][i,]))
          pvals = c(pvals, p)
        }
      }
      # store the pvals in pval_groups list
      names(pvals) = rownames(realenrich[[c]])
      pval_groups[[c]] = pvals
    }
    return(pval_groups)
  }
}


