library(RERconverge)
rerpath = find.package('RERconverge')

#read trees
toytreefile = "subsetMammalGeneTrees.txt"
toyTrees=readTrees(paste(rerpath,"/extdata/",toytreefile,sep=""), max.read = 200)

#### Work on permulation for specific trees with internal foregrounds with no tip foregrounds
fg_vec_test2 = c("Killer_whale", "Dolphin", "Walrus", "Seal", "Manatee")
sisters_list_test2 = list("clade1"=c("Killer_whale", "Dolphin"))
root_sp = "Human"

fgTree = foreground2TreeClades(fg_vec_test2,sisters_list_test2,toyTrees,plotTree=T)
# zero out Killer_whale and Dolphin
rm_fgs = c('Killer_whale', 'Dolphin')
idx_rm_fgs = which(fgTree$tip.label %in% rm_fgs)
idx_fg_edges = which(fgTree$edge[,2] %in% idx_rm_fgs)
fgTree$edge.length[idx_fg_edges] = 0
plot(fgTree)


## Get independent lineage structure in fgTree
idx_fg_edges = which(fgTree$edge.length ==1)
fg_edges = fgTree$edge[idx_fg_edges,]

all_edges = fgTree$edge

unq_start_nodes = unique(fg_edges[,1])
unq_end_nodes = unique(fg_edges[,2])

# remove nested starting points
nested_nodes = intersect(unq_start_nodes, unq_end_nodes)
if (length(nested_nodes) > 1){
  unq_start_nodes = unq_start_nodes[-(which(unq_start_nodes %in% nested_nodes))]
  unq_end_nodes = unq_end_nodes[-(which(unq_end_nodes %in% nested_nodes))]
}

# delineate the depth structure of each lineage
lineages=list()
counter = 0
for (i in 1:length(unq_start_nodes)){
  start_node_i = unq_start_nodes[i]
  # identify foreground branches that start from this node
  idx_fg_branches_i = which(fg_edges[,1] == start_node_i)
  target_nodes_i = fg_edges[idx_fg_branches_i, 2]
  for (j in 1:length(target_nodes_i)){
    target_node_j = target_nodes_i[j]
    if (target_node_j <= length(fgTree$tip.label)){
      counter = counter+1
      lineage_j = data.frame("depth_order"= 0, "phenval"=1)
      lineages[[counter]] = lineage_j
    } else {
      tip_node = FALSE
      counter = counter+1
      depth_order = 1
      while (tip_node){


      }
    }
  }
}





getDepthOrder(fgTree)

fg_tree = fgTree

getForegroundInfoFromTree=function(fg_tree,trees,plotTree=T,useSpecies=NULL){
  if (length(useSpecies)==0){
    useSpecies = trees$masterTree$tip.label
  }

  # Start with a temp phenotype tree assuming that all internal nodes are foregrounds
#  fg_tree = foreground2Tree(fg_vec,trees,plotTree=F,clade="all",useSpecies=useSpecies)


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
      for (nn in 1:(length(nodeId.ca)-1)){
        if (nodeId.ca[nn] == nodeId.ca[nn+1]){
          nodeId.desc = nodeIds.fg.edge[which(nodeIds.fg.edge[,1]==nodeId.ca[nn]),2]
          if (length(which(nodeId.desc %in% tip.sisters)) > 0){
            fg_ca = c(fg_ca,nodeId.ca[nn])
            fg.sisters.table = rbind(fg.sisters.table, nodeIds.fg.edge[which(nodeIds.fg.edge[,1]==nodeId.ca[nn]),2])
            nodes_addressed = c(nodes_addressed, nodeId.ca[nn])
          } else {
            if (length(which(trees$masterTree$tip.label[nodeId.desc] %in% fg_vec)) == 2){
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
    rownames(fg.sisters.table) = fg_ca
  }

  if (plotTree==T){
    plot(fg_tree)
  }
  output.list = list("fg.sisters.table"=fg.sisters.table,"tree"=fg_tree)
  output.list
}



### test foreground2TreeClades
fg_vec_test1 = c("Killer_whale", "Dolphin", "Walrus", "Seal", "Manatee", 'Sheep','Goat','Tibetan_antelope', 'Cow', 'Pig',
           'Rhinoceros', 'Horse')
sisters_list_test1 = list("clade1"=c("Killer_whale", "Dolphin"),
                    'clade4'=c('Sheep', 'Goat'), 'clade5'=c('Tibetan_antelope', 'clade4'), 'clade6'=c('clade5', 'Cow'),
                    'clade7'=c('clade6', 'clade1'),
                    'clade8'=c('clade7','Pig'),
                    'clade9'=c('Rhinoceros', 'Horse'))
root_sp = "Human"

fgTree_test1 = foreground2TreeClades(fg_vec_test1,sisters_list_test1,toyTrees,plotTree=F)
fgplot1 = plotTreeHighlightBranches(fgTree_test1,
                                        hlspecies=which(fgTree_test1$edge.length==1),
                                        hlcols="blue", main="Marine mammals trait tree")
out_test = getBinaryPermulationInputsFromTree(fgTree_test1)
test_tree = foreground2TreeClades(out_test$fg_vec, out_test$sisters_list, toyTrees, plotTree=T)



### test binary permulation functions
fg_vec_test2 = c("Killer_whale", "Dolphin", "Walrus", "Seal", "Manatee")
sisters_list_test2 = list("clade1"=c("Killer_whale", "Dolphin"))
root_sp = "Human"

fgTree_test2 = foreground2TreeClades(fg_vec_test2,sisters_list_test2,toyTrees,plotTree=T)
pathvec_test2 = tree2PathsClades(fgTree_test2, toyTrees)

test_bin1 = simBinPhenoCC(toyTrees, toyTrees$masterTree, root_sp, fg_vec_test2, sisters_list_test2, pathvec_test2, plotTreeBool=T)


fg_vec_test3 = c("Killer_whale", "Dolphin", "Walrus", "Seal", "Manatee", "Ferret")
sisters_list_test3 = list("clade1"=c("Killer_whale", "Dolphin"),
                          "clade2"=c("Walrus", "Seal"),
                          "clade3"=c("clade2", "Ferret"))
root_sp = "Human"

fgTree_test3 = foreground2TreeClades(fg_vec_test3,sisters_list_test3,toyTrees,plotTree=T)
pathvec_test3 = tree2PathsClades(fgTree_test3, toyTrees)

test_bin3 = simBinPhenoCC(toyTrees, toyTrees$masterTree, root_sp, fg_vec_test3, sisters_list_test3, pathvec_test3, plotTreeBool=T)

plot(fgTree_test3)



fg_vec_test4 = c("Killer_whale", "Dolphin", "Walrus", "Seal", "Manatee", "Ferret", "Panda")
sisters_list_test4 = list("clade1"=c("Killer_whale", "Dolphin"),
                          "clade2"=c("Walrus", "Seal"),
                          "clade3"=c("clade2", "Ferret"),
                          "clade4"=c("clade3", "Panda"))
root_sp = "Human"

fgTree_test4 = foreground2TreeClades(fg_vec_test4,sisters_list_test4,toyTrees,plotTree=T)
plot(fgTree_test4)

pathvec_test4 = tree2PathsClades(fgTree_test4, toyTrees)

test_bin4 = simBinPhenoCC(toyTrees, toyTrees$masterTree, root_sp, fg_vec_test4, sisters_list_test4, pathvec_test4, plotTreeBool=T)




fg_vec_test5 = c("Killer_whale", "Dolphin", "Walrus", "Seal", "Manatee", "Ferret", "Panda")
sisters_list_test5 = NULL
root_sp = "Human"

fgTree_test5 = foreground2Tree(fg_vec_test5,toyTrees,plotTree=F, clade='terminal')
plot(fgTree_test5)

pathvec_test5 = tree2PathsClades(fgTree_test5, toyTrees)

test_bin5 = simBinPhenoCC(toyTrees, toyTrees$masterTree, root_sp, fg_vec_test5, sisters_list_test5, pathvec_test5, plotTreeBool=T)



#### debugging simBinPhenoCC with getBinaryPermulationInputsFromTree
trees = toyTrees
mastertree=trees$masterTree
fg_vec = fg_vec_test5
sisters_list = sisters_list_test5
pathvec = pathvec_test5

test_out_4 = generatePermulatedBinPhen(trees$mastertree, 3, trees,
                                       root_sp, fg_vec_test4, sisters_list_test4,
                                       pathvec_test4, permmode="cc")

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
rep_tree = function(num_input,tree){
  return(tree)
}



simBinPhenoCC=function(trees, mastertree, root_sp, fg_vec, sisters_list=NULL, pathvec, plotTreeBool=F){
  tip.labels = mastertree$tip.label
  res = getForegroundInfoClades(fg_vec,sisters_list,trees,plotTree=F,useSpecies=tip.labels)  #### This function has problems
  fg_tree = res$tree
  fg.table = res$fg.sisters.table

  if (!is.null(sisters_list)){
    fg_tree_info = getBinaryPermulationInputsFromTree(fg_tree)
    num_tip_sisters_true = unlist(fg_tree_info$sisters_list)
    num_tip_sisters_true = num_tip_sisters_true[which(num_tip_sisters_true %in% tip.labels)]
    num_tip_sisters_true = length(num_tip_sisters_true)
    fg_tree_depth_order = getDepthOrder(fg_tree)
  }
  t=root.phylo(trees$masterTree, root_sp, resolve.root = T)
  rm=ratematrix(t, pathvec)

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
      t=foreground2Tree(top, trees, clade="all", plotTree = F)
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
  if (plotTreeBool){
    plot(t)
  }
  return(t)
}






