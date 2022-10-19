library(RERconverge)
rerpath = find.package('RERconverge')

#read trees
toytreefile = "subsetMammalGeneTrees.txt"
trees=readTrees(paste(rerpath,"/extdata/",toytreefile,sep=""), max.read = 200)

fg_vec = c("Killer_whale", "Dolphin", "Walrus", "Seal", "Manatee", 'Sheep','Goat','Tibetan_antelope', 'Cow', 'Pig',
           'Rhinoceros', 'Horse')
sisters_list = list("clade1"=c("Killer_whale", "Dolphin"),
                    'clade4'=c('Sheep', 'Goat'), 'clade5'=c('Tibetan_antelope', 'clade4'), 'clade6'=c('clade5', 'Cow'),
                    'clade7'=c('clade6', 'clade1'),
                    'clade8'=c('clade7','Pig'),
                    'clade9'=c('Rhinoceros', 'Horse'))
root_sp = "Human"


#trees = toyTrees
mastertree = trees$masterTree


fgTree = foreground2TreeClades(fg_vec,sisters_list,toyTrees,plotTree=F)
fgplot1 = plotTreeHighlightBranches(fgTree,
                                        hlspecies=which(fgTree$edge.length==1),
                                        hlcols="blue", main="Marine mammals trait tree")


out_test = getPhenotypePermulationInputsFromTree(fgTree)

test_tree = foreground2TreeClades(out_test$fg_vec, out_test$sisters_list, trees, plotTree=T)

getPhenotypePermulationInputsFromTree=function(fgTree){
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
            node_daughter = daughters_j[which(daughters_j > length(mastertree$tip.label))]
            if (node_daughter %in% order_assigned){
              depth_order[idx_na_j] = 2
              order_assigned = c(order_assigned, parent_j)
            }
          } else if (num_tip_daughters==0){
            node_daughters = daughters_j[which(daughters_j > length(mastertree$tip.label))]
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
        if (length(which(daughters_i <= length(mastertree$tip.label))) == 2){
          counter = counter+1
          tip_daughters = mastertree$tip.label[daughters_i]
          sisters_list[[counter]] = tip_daughters
          names(sisters_list)[counter] = names(daughters_info_order_j)[i]
          nodes_addressed = c(nodes_addressed, names(daughters_info_order_j)[i])
        } else if (length(which(daughters_i <=length(mastertree$tip.label))) == 1){
          counter = counter+1
          tip_daughter_id = daughters_i[which(daughters_i <= length(mastertree$tip.label))]
          tip_daughter = mastertree$tip.label[tip_daughter_id]
          node_daughter_id = daughters_i[which(daughters_i > length(mastertree$tip.label))]
          sisters_list[[counter]] = c(node_daughter_id, tip_daughter)
          names(sisters_list)[counter] = names(daughters_info_order_j)[i]
          nodes_addressed = c(nodes_addressed, names(daughters_info_order_j)[i])
        } else if (length(which(daughters_i <=length(mastertree$tip.label))) == 0){
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



### get fg_vec and sisters_list from fgTree
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
          node_daughter = daughters_j[which(daughters_j > length(mastertree$tip.label))]
          if (node_daughter %in% order_assigned){
            depth_order[idx_na_j] = 2
            order_assigned = c(order_assigned, parent_j)
          }
        } else if (num_tip_daughters==0){
          node_daughters = daughters_j[which(daughters_j > length(mastertree$tip.label))]
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
      if (length(which(daughters_i <= length(mastertree$tip.label))) == 2){
        counter = counter+1
        tip_daughters = mastertree$tip.label[daughters_i]
        sisters_list[[counter]] = tip_daughters
        names(sisters_list)[counter] = names(daughters_info_order_j)[i]
        nodes_addressed = c(nodes_addressed, names(daughters_info_order_j)[i])
      } else if (length(which(daughters_i <=length(mastertree$tip.label))) == 1){
        counter = counter+1
        tip_daughter_id = daughters_i[which(daughters_i <= length(mastertree$tip.label))]
        tip_daughter = mastertree$tip.label[tip_daughter_id]
        node_daughter_id = daughters_i[which(daughters_i > length(mastertree$tip.label))]
        sisters_list[[counter]] = c(node_daughter_id, tip_daughter)
        names(sisters_list)[counter] = names(daughters_info_order_j)[i]
        nodes_addressed = c(nodes_addressed, names(daughters_info_order_j)[i])
      } else if (length(which(daughters_i <=length(mastertree$tip.label))) == 0){
        counter = counter+1
        sisters_list[[counter]] = daughters_i
        names(sisters_list)[counter] = names(daughters_info_order_j)[i]
        nodes_addressed = c(nodes_addressed, names(daughters_info_order_j)[i])
      }
    }

  }

}


fgOut = foreground2TreeClades(tip_foregrounds,sisters_list,toyTrees,plotTree=F)
fgplotOut = plotTreeHighlightBranches(fgOut,
                                    hlspecies=which(fgOut$edge.length==1),
                                    hlcols="blue", main="Marine mammals trait tree")

#### Debugging foreground2TreeClades


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
    all.nodeId.ca = sort(nodeIds.fg.edge[,1])
    count_all_nodeId_ca = table(all.nodeId.ca)
    unq.nodeId.ca = unique(all.nodeId.ca)
    fg_ca = vector("integer",length=0) # node IDs of the common ancestor foregrounds
    nodes_addressed = NULL
    counter=0
    while (length(unq.nodeId.ca) != length(nodes_addressed)){
      counter = counter+1
      print(counter)
      nodeId.ca = sort(all.nodeId.ca[which(!(all.nodeId.ca %in% nodes_addressed))])
      for (nn in 1:(length(nodeId.ca)-1)){
        print(paste('Evaluating', nodeId.ca[nn]))
        if (nodeId.ca[nn] == nodeId.ca[nn+1]){
          print('- node has 2 foreground daughters in tmp tree')
          nodeId.desc = nodeIds.fg.edge[which(nodeIds.fg.edge[,1]==nodeId.ca[nn]),2]
          if (length(which(nodeId.desc %in% tip.sisters)) > 0){
            print('--- node has 1 or 2 tip foregrounds in real tree, node is an ancestral foreground')
            fg_ca = c(fg_ca,nodeId.ca[nn])
            fg.sisters.table = rbind(fg.sisters.table, nodeIds.fg.edge[which(nodeIds.fg.edge[,1]==nodeId.ca[nn]),2])
            nodes_addressed = c(nodes_addressed, nodeId.ca[nn])
          } else {
            print('--- node has no tip foreground daughters in real tree; node either has 2 ancestral foreground daughters or node is not an ancestral foreground')
            if (length(which(mastertree$tip.label[nodeId.desc] %in% fg_vec)) == 2){
              print('----- node is not an ancestral foreground')
              fg_tree$edge.length[which(edge[,2]==nodeId.ca[nn])] = 0
              nodes_addressed = c(nodes_addressed, nodeId.ca[nn])
            } else {
              if (length(which(nodeId.desc %in% nodes_addressed)) == 2){
                print('----- node is an ancestral foreground with 2 ancestral daughters')
                fg_ca = c(fg_ca,nodeId.ca[nn])
                fg.sisters.table = rbind(fg.sisters.table, nodeIds.fg.edge[which(nodeIds.fg.edge[,1]==nodeId.ca[nn]),2])
                nodes_addressed = c(nodes_addressed, nodeId.ca[nn])
              }
            }
          }
        } else {
          print('- node has 1 foreground daughter in tmp tree')
          nodeId.desc = nodeIds.fg.edge[which(nodeIds.fg.edge[,1]==nodeId.ca[nn]),2]
          if (length(nodeId.desc) == 2){
            if (nodeId.ca[nn] != nodeId.ca[nn-1]){
              fg_tree$edge.length[which(edge[,2] == nodeId.ca[nn])] = 0
              nodes_addressed = c(nodes_addressed, nodeId.ca[nn])
              nodes_addressed = unique(nodes_addressed)
            }
          } else {
            print('--- node is a parent of 1 independent branch in tmp tree')
            nodes_addressed = c(nodes_addressed, nodeId.ca[nn])
          }
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

foreground2TreeClades=function(fg_vec,sisters_list=NULL,trees,plotTree=T,useSpecies=NULL){
  res.list = getForegroundInfoClades(fg_vec,sisters_list=sisters_list,trees,plotTree=plotTree,useSpecies=useSpecies)
  fg.tree = res.list$tree
  fg.tree
}




####
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



pathvec = tree2PathsClades(fg_tree, trees)



simBinPhenoRank=function(trees, root_sp, fg_vec, sisters_list=NULL, plotTreeBool=F){
  mastertree = trees$masterTree
  tip.labels = mastertree$tip.label
  res = getForegroundInfoClades(fg_vec,sisters_list,trees,plotTree=F,useSpecies=tip.labels)
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


treetest = simBinPhenoRank(trees, root_sp, fg_vec, sisters_list, plotTreeBool=T)
plot(treetest)




sims=sim.char(t, rm, nsim = 1)
nam=rownames(sims)
s=as.data.frame(sims)
simulatedvec=s[,1]
names(simulatedvec)=nam
top=names(sort(simulatedvec, decreasing = TRUE))[1:tips]
t=foreground2Tree(top, trees, clade="all", plotTree = F)
blsum=sum(t$edge.length)








geo <- get(data(geospiza))

## Continuous character -- univariate
usims <- sim.char(geo$phy, 0.02, 100)

## Use a simulated dataset in fitContinuous()
fitC <- fitContinuous(geo$phy, usims[,,1], model="BM", control=list(niter=10), ncores=2)

## Continuous character -- multivariate
s <- ratematrix(geo$phy, geo$dat)
csims <- sim.char(geo$phy, s, 100)

## Discrete character -- univariate
q <- list(rbind(c(-.5, .5), c(.5, -.5)))
dsims <- sim.char(geo$phy, q, model="discrete", n=10)

## Use a simulated dataset in fitDiscrete()
fitD <- fitDiscrete(geo$phy, dsims[,,1], model="ER", niter=10, ncores=2)

## Discrete character -- multivariate
qq <- list(rbind(c(-.5, .5), c(.5, -.5)), rbind(c(-.05, .05), c(.05, -.05)))
msims <- sim.char(geo$phy, qq, model="discrete", n=10)




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




