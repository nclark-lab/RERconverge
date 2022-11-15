library(RERconverge)
rerpath = find.package('RERconverge')

#read trees
toytreefile = "subsetMammalGeneTrees.txt"
toyTrees=readTrees(paste(rerpath,"/extdata/",toytreefile,sep=""), max.read = 200)
#test_read=readLines(paste(rerpath,"/extdata/",toytreefile,sep=""))

#### Work on permulation for specific trees with internal foregrounds with no tip foregrounds
fg_vec_test2 = c("Killer_whale", "Dolphin", "Walrus", "Seal", "Manatee")
sisters_list_test2 = list("clade1"=c("Killer_whale", "Dolphin"))
root_sp = "Human"

#fgTree = foreground2TreeClades(fg_vec_test2,sisters_list_test2,toyTrees,plotTree=T)

# drop 2 species to test SSM
#trees = toyTrees
tree = toyTrees$masterTree
tree = drop.tip(tree, c('Ferret', 'Killer_whale'))


fgTree_test2 = foreground2TreeClades(fg_vec_test2,sisters_list_test2,toyTrees,plotTree=T)
pathvec_test2 = tree2PathsClades(fgTree_test2, toyTrees)

test_bin2 = simBinPhenoCC(toyTrees, toyTrees$masterTree, root_sp, fg_vec_test2, sisters_list_test2, pathvec_test2, plotTreeBool=T)
test_ssm2 = simBinPhenoSSM(tree, toyTrees, root_sp, fg_vec_test2, sisters_list=sisters_list_test2, pathvec_test2, plotTreeBool=T)

idx_tree = 85
tree = toyTrees$trees[[idx_tree]]
print(length(tree$tip.label))
print(length(which(fg_vec_test2 %in% tree$tip.label)))
simBinPhenoSSM(toyTrees$trees[[3]], toyTrees, root_sp, fg_vec_test2, sisters_list=sisters_list_test2, pathvec_test2, plotTreeBool=T)

### test generatePermulatedBinPhen
test_permulated_trees_cc = generatePermulatedBinPhen(toyTrees$masterTree, 10, toyTrees, root_sp, fg_vec_test2, sisters_list_test2, pathvec_test2, 'cc')
test_permulated_trees_ssm = generatePermulatedBinPhen(tree, 10, toyTrees, root_sp, fg_vec_test2, sisters_list_test2, pathvec_test2, 'ssm')

### test generatePermulatedBinPhenSSMBatched
test_permulated_trees_ssm_batched = generatePermulatedBinPhenSSMBatched(toyTrees$trees, 10, toyTrees, root_sp, fg_vec_test2, sisters_list_test2, pathvec_test2)

trees_list = toyTrees$trees
numperms=5
trees = toyTrees
fg_vec = fg_vec_test2
sisters_list = sisters_list_test2
pathvec = pathvec_test2

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



