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

### test generatePermulatedBinPhen ####
test_permulated_trees_cc = generatePermulatedBinPhen(toyTrees$masterTree, 10, toyTrees, root_sp, fg_vec_test2, sisters_list_test2, pathvec_test2, 'cc')
test_permulated_trees_ssm = generatePermulatedBinPhen(tree, 10, toyTrees, root_sp, fg_vec_test2, sisters_list_test2, pathvec_test2, 'ssm')

### test generatePermulatedBinPhenSSMBatched ####
test_permulated_trees_ssm_batched = generatePermulatedBinPhenSSMBatched(toyTrees$trees, 10, toyTrees, root_sp, fg_vec_test2, sisters_list_test2, pathvec_test2)

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



### test getPermsBinary
numperms = 5
fg_vec = fg_vec_test2
sisters_list = sisters_list_test2
trees = toyTrees
mastertree = trees$masterTree
permmode='cc'
method='k'
min.pos=2
trees_list=NULL
calculateenrich=F
annotlist=NULL

RERmat = getAllResiduals(toyTrees,transform = "sqrt", weighted = T, scale = T)

test_get_perms_binary_cc = getPermsBinary(10, fg_vec_test2, sisters_list_test2, root_sp, RERmat,
                                          toyTrees, toyTrees$mastertree, permmode="cc", method="k",
                                          min.pos=2, trees_list=NULL, calculateenrich=F, annotlist=NULL)
test_get_perms_binary_ssm = getPermsBinary(10, fg_vec_test2, sisters_list_test2, root_sp, RERmat,
                                          toyTrees, toyTrees$mastertree, permmode="ssm", method="k",
                                          min.pos=2, trees_list=NULL, calculateenrich=F, annotlist=NULL)

realcor=correlateWithBinaryPhenotype(RERmat, pathvec_test2)

### test permpvalcor
permvals = test_get_perms_binary_cc

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
    permstats[count] = -log10(permpvals[count])*(realstat[count]-median_permcor)
  }
  count=count+1
}




######## getPermsBinary SSM
numperms=10
fg_vec = fg_vec_test2
sisters_list = sisters_list_test2
trees = toyTrees
mastertree=trees$masterTree
permmode='ssm'
method='k'
min.pos=2
trees_list = NULL

if (is.null(trees_list)){
  trees_list = trees$trees
}
RERmat_true = RERmat


pathvec = foreground2Paths(fg_vec, trees, clade="all",plotTree=F)
col_labels = colnames(trees$paths)
names(pathvec) = col_labels


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



