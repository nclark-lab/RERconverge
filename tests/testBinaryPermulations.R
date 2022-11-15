library(RERconverge)
rerpath = find.package('RERconverge')

#read trees
toytreefile = "subsetMammalGeneTrees.txt"
toyTrees=readTrees(paste(rerpath,"/inst/extdata/",toytreefile,sep=""), max.read = 200)
#test_read=readLines(paste(rerpath,"/extdata/",toytreefile,sep=""))

#### Work on permulation for specific trees with internal foregrounds with no tip foregrounds
fg_vec_test2 = c("Killer_whale", "Dolphin", "Walrus", "Seal", "Manatee")
sisters_list_test2 = list("clade1"=c("Killer_whale", "Dolphin"))
root_sp = "Human"

fgTree = foreground2TreeClades(fg_vec_test2,sisters_list_test2,toyTrees,plotTree=T)

# drop 2 species to test SSM
trees = toyTrees
tree = trees$masterTree
tree = drop.tip(tree, c('Ferret', 'Killer_whale'))


fgTree_test2 = foreground2TreeClades(fg_vec_test2,sisters_list_test2,toyTrees,plotTree=T)
pathvec_test2 = tree2PathsClades(fgTree_test2, toyTrees)


simBinPhenoSSM(tree, toyTrees, root_sp, fg_vec_test2, sisters_list=sisters_list_test2, pathvec_test2, plotTreeBool=T)
simBinPhenoSSM(toyTrees$masterTree, toyTrees, root_sp, fg_vec_test2, sisters_list=sisters_list_test2, pathvec_test2, plotTreeBool=T)

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
    print('here')
#    stop('stopped')
    fg_k = tip.labels[ind_fg] # the list of the observed foreground animals that exist in the gene tree

    res = getForegroundInfoClades(fg_k,sisters_list,trees,plotTree=F,useSpecies=tip.labels)
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
        t=foreground2Tree(top, trees, clade="all", plotTree = F, useSpecies=tip.labels)
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

fg_vec = fg_k
useSpecies=tip.labels

getForegroundInfoClades=function(fg_vec,sisters_list=NULL,trees,plotTree=T,useSpecies=NULL){
  if (length(useSpecies)==0){
    useSpecies = trees$masterTree$tip.label
  }

  if (is.null(sisters_list)){
    fg.sisters.table=NULL
    fg_tree = foreground2Tree(fg_vec, trees, plotTree=F, clade="terminal", useSpecies=useSpecies)
  } else {
    # Start with a temp phenotype tree assuming that all internal nodes are foregrounds
    fg_tree = foreground2Tree(fg_vec,trees,plotTree=F,clade="all",useSpecies=useSpecies)
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


fgTree = fg_tree


getDepthOrder=function(fgTree){
  unq_edge_lengths = unique(fgTree$edge.length)
  if (length(which(!(unq_edge_lengths %in% c(0,1)))) > 0){
    stop('Phenotype must be binary.')
  }

  idx_fg_branches = which(fgTree$edge.length == 1)
  fg_edges = fgTree$edge[idx_fg_branches,]
  all_edges = fgTree$edge

  num_tip_species = length(fgTree$tip.label)
  tip_fg_edges = fg_edges[which(fg_edges[,2] <= num_tip_species),]
  tip_foregrounds = fgTree$tip.label[tip_fg_edges[,2]]
  node_fg_edges = fg_edges[which(fg_edges[,2] > num_tip_species),]

  idx_node_edges = which(fg_edges[,2] > num_tip_species)
  if (length(idx_node_edges) == 1){
    node_fg_edges = fg_edges[which(fg_edges[,2] > num_tip_species),]
    node_fg_edges = t(as.data.frame(node_fg_edges))
  }
  if (length(idx_node_edges) == 0) {
    sisters_list = NULL
    depth_order=NULL
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




getBinaryPermulationInputsFromTree=function(fgTree){
  unq_edge_lengths = unique(fgTree$edge.length)
  if (length(which(!(unq_edge_lengths %in% c(0,1)))) > 0){
    stop('Phenotype must be binary.')
  }

  idx_fg_branches = which(fgTree$edge.length == 1)
  fg_edges = fgTree$edge[idx_fg_branches,]
  all_edges = fgTree$edge

  num_tip_species = length(fgTree$tip.label)
  tip_fg_edges = fg_edges[which(fg_edges[,2] <= num_tip_species),]
  tip_foregrounds = fgTree$tip.label[tip_fg_edges[,2]]

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

tree = drop.tip(toyTrees$masterTree, c('Walrus', 'Ferret', 'Dolphin', 'Killer_whale'))

test_bin1 = simBinPhenoCC(toyTrees, toyTrees$masterTree, root_sp, fg_vec_test2, sisters_list_test2, pathvec_test2, plotTreeBool=T)
test_ssm1 = simBinPhenoSSM(tree, toyTrees, root_sp, fg_vec_test2, sisters_list=sisters_list_test2, pathvec_test2, plotTreeBool=T)



fg_vec_test3 = c("Killer_whale", "Dolphin", "Walrus", "Seal", "Manatee", "Ferret")
sisters_list_test3 = list("clade1"=c("Killer_whale", "Dolphin"),
                          "clade2"=c("Walrus", "Seal"),
                          "clade3"=c("clade2", "Ferret"))
root_sp = "Human"

fgTree_test3 = foreground2TreeClades(fg_vec_test3,sisters_list_test3,toyTrees,plotTree=T)
pathvec_test3 = tree2PathsClades(fgTree_test3, toyTrees)
tree = drop.tip(toyTrees$masterTree, c('Walrus', 'Ferret', 'Seal'))

test_bin3 = simBinPhenoCC(toyTrees, toyTrees$masterTree, root_sp, fg_vec_test3, sisters_list_test3, pathvec_test3, plotTreeBool=T)
test_ssm3 = simBinPhenoSSM(tree, toyTrees, root_sp, fg_vec_test3, sisters_list=sisters_list_test3, pathvec_test3, plotTreeBool=T)




fg_vec_test4 = c("Killer_whale", "Dolphin", "Walrus", "Seal", "Manatee", "Ferret", "Panda")
sisters_list_test4 = list("clade1"=c("Killer_whale", "Dolphin"),
                          "clade2"=c("Walrus", "Seal"),
                          "clade3"=c("clade2", "Ferret"),
                          "clade4"=c("clade3", "Panda"))
root_sp = "Human"

fgTree_test4 = foreground2TreeClades(fg_vec_test4,sisters_list_test4,toyTrees,plotTree=T)
pathvec_test4 = tree2PathsClades(fgTree_test4, toyTrees)
tree = drop.tip(toyTrees$masterTree, c('Panda', 'Walrus', 'Seal'))

test_bin4 = simBinPhenoCC(toyTrees, toyTrees$masterTree, root_sp, fg_vec_test4, sisters_list_test4, pathvec_test4, plotTreeBool=T)
test_ssm4 = simBinPhenoSSM(tree, toyTrees, root_sp, fg_vec_test4, sisters_list=sisters_list_test4, pathvec_test4, plotTreeBool=T)




fg_vec_test5 = c("Killer_whale", "Dolphin", "Walrus", "Seal", "Manatee", "Ferret", "Panda")
sisters_list_test5 = NULL
root_sp = "Human"

fgTree_test5 = foreground2Tree(fg_vec_test5,toyTrees,plotTree=F, clade='terminal')
plot(fgTree_test5)

pathvec_test5 = tree2PathsClades(fgTree_test5, toyTrees)

test_bin5 = simBinPhenoCC(toyTrees, toyTrees$masterTree, root_sp, fg_vec_test5, sisters_list_test5, pathvec_test5, plotTreeBool=T)



