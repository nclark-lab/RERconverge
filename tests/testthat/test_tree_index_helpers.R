test_that("prepareTreeForTT aligns subset mammal trees with master", {
  subset <- mammal_subset_trees(2)
  master <- TreeTools::Preorder(TreeTools::SortTree(subset[[1]]))
  gene <- subset[[2]]
  prepared <- RERconverge:::prepareTreeForTT(gene, master)
  expect_true(ape::is.rooted(prepared))
  expect_setequal(prepared$tip.label, intersect(master$tip.label, gene$tip.label))

  map <- RERconverge:::matchAllnodesTT(prepared, master)
  expect_equal(nrow(map), length(prepared$tip.label) + prepared$Nnode)
  expect_true(all(map[, 2] %in% seq_len(master$Nnode + length(master$tip.label))))
})

test_that("allPathsTT matches ancestor distances and matrix indices", {
  tree <- ape::read.tree(text = "((A:1,B:2):3,(C:4,D:5):6);")
  tree <- TreeTools::Preorder(TreeTools::SortTree(tree))
  paths <- RERconverge:::allPathsTT(tree)
  dmat <- ape::dist.nodes(tree)
  idx_vals <- vapply(seq_len(nrow(paths$nodeId)), function(i) {
    dmat[paths$nodeId[i, 1], paths$nodeId[i, 2]]
  }, numeric(1))
  expect_equal(paths$dist, idx_vals)

  nn <- length(tree$tip.label) + tree$Nnode
  has_path <- matrix(FALSE, nn, nn)
  for (i in seq_len(nn)) {
    anc <- RERconverge:::getAncestors(tree, i)
    if (length(anc) > 0) {
      has_path[i, anc] <- TRUE
    }
  }
  expect_equal(sort(paths$matIndex[has_path]), seq_along(paths$dist))
  expect_true(all(is.na(paths$matIndex[!has_path])))
})

test_that("namePathsWSpeciesTT labels tip-directed columns", {
  trees <- build_simple_trees()
  cnames <- RERconverge:::namePathsWSpeciesTT(trees)
  expect_setequal(na.omit(cnames[cnames != ""]), trees$masterTree$tip.label)
})

test_that("allPathsMasterRelativeTT respects pruning and detects discordance", {
  trees <- build_simple_trees()
  master <- trees$masterTree
  masterPaths <- RERconverge:::allPathsTT(master)
  full_vals <- RERconverge:::allPathsMasterRelativeTT(master, master, masterPaths)
  expect_equal(full_vals, masterPaths$dist)

  discordant <- ape::read.tree(text = "((A:1,C:1):2,(B:1,D:1):2);")
  res <- suppressWarnings(RERconverge:::allPathsMasterRelativeTT(discordant, master, masterPaths))
  expect_true(all(is.na(res)))
})
