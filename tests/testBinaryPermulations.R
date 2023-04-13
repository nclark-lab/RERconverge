library(RERconverge)
rerpath = find.package('RERconverge')

#read trees
toytreefile = "subsetMammalGeneTrees.txt"
toyTrees=readTrees(paste(rerpath,"/extdata/",toytreefile,sep=""), max.read = 200)
#test_read=readLines(paste(rerpath,"/extdata/",toytreefile,sep=""))

### test midpoint rooting
simBinPhenoCCmidpoint=function(trees, mastertree, fg_vec, sisters_list=NULL, pathvec, plotTreeBool=F){
  tip.labels = mastertree$tip.label
  res = getForegroundInfoClades(fg_vec,sisters_list,trees,plotTree=F,useSpecies=tip.labels)  #### This function has problems
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
      t_iter=foreground2Tree(top, trees, clade="all", plotTree = F)
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



### no internal foregrounds
fg_vec = c("Killer_whale", "Dolphin", "Walrus", "Seal", "Manatee")
sisters_list = NULL

fgTree = foreground2TreeClades(fg_vec,sisters_list,toyTrees,plotTree=T)
pathvec = tree2PathsClades(fgTree, toyTrees)
permphen = simBinPhenoCCmidpoint(toyTrees, toyTrees$masterTree, fg_vec, sisters_list, pathvec, plotTreeBool=T)

### with internal foreground
sisters_list_v2 = list("clade1"=c("Killer_whale", "Dolphin"))
fgTree = foreground2TreeClades(fg_vec,sisters_list_v2,toyTrees,plotTree=T)
pathvec = tree2PathsClades(fgTree, toyTrees)
permphen = simBinPhenoCCmidpoint(toyTrees, toyTrees$masterTree, fg_vec, sisters_list_v2, pathvec, plotTreeBool=T)




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

### test generatePermulatedBinPhen ####
test_permulated_trees_cc = generatePermulatedBinPhen(toyTrees$masterTree, 10, toyTrees, root_sp, fg_vec_test2, sisters_list_test2, pathvec_test2, 'cc')
test_permulated_trees_ssm = generatePermulatedBinPhen(tree, 10, toyTrees, root_sp, fg_vec_test2, sisters_list_test2, pathvec_test2, 'ssm')

### test generatePermulatedBinPhenSSMBatched ####
test_permulated_trees_ssm_batched = generatePermulatedBinPhenSSMBatched(toyTrees$trees, 10, toyTrees, root_sp, fg_vec_test2, sisters_list_test2, pathvec_test2)

test_permulated_trees_ssm_batched=permulated.binphens
# check if all the resulting permulated trees have the same number of foregrounds as in the observed
for (i in 1:length(toyTrees$trees)){
  print(i)
  tree_i = toyTrees$trees[[1]]
  tips_observed = tree_i$tip.label
  fg_observed = length(which(tips_observed %in% fg_vec))
  ssm_trees = test_permulated_trees_ssm_batched[[1]]
  for (j in 1:length(ssm_trees)){
    ssm_tree_j = ssm_trees[[j]]
    tips_ssm_j = ssm_tree_j$tip.label
    if (length(tips_observed)!=length(tips_ssm_j)){
      stop('mismatch in number of tip species')
    }
    idx_fg = which(ssm_tree_j$edge.length == 1)
    fg_edges = ssm_tree_j$edge[idx_fg,]
    idx_tip_fgs = which(fg_edges[,2] <= length(tips_ssm_j))
    if (length(idx_tip_fgs) != fg_observed){
      stop('mismatch in number of tip foregrounds')
    }
  }
}
# all the generated SSM trees correctly match the corresponding gene trees



RERmat = getAllResiduals(toyTrees,transform = "sqrt", weighted = T, scale = T)

test_get_perms_binary_cc = getPermsBinary(10, fg_vec_test2, sisters_list_test2, root_sp, RERmat,
                                          toyTrees, toyTrees$mastertree, permmode="cc", method="k",
                                          min.pos=2, trees_list=NULL, calculateenrich=F, annotlist=NULL)
test_get_perms_binary_ssm = getPermsBinary(100, fg_vec_test2, sisters_list_test2, root_sp, RERmat,
                                          toyTrees, toyTrees$mastertree, permmode="ssm", method="k",
                                          min.pos=2, trees_list=NULL, calculateenrich=F, annotlist=NULL)

realcor=correlateWithBinaryPhenotype(RERmat, pathvec_test2)

test_ssm_pval = permpvalcor(realcor, test_get_perms_binary_ssm)
test_cc_pval = permpvalcor(realcor, test_get_perms_binary_cc)









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



