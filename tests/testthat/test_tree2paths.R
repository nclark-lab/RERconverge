test_that("tree2Paths reproduces TT path ordering and honors binarize/categorical options", {
  trees <- build_simple_trees()
  master <- trees$masterTree
  master <- TreeTools::Preorder(TreeTools::SortTree(master))
  pheno <- master
  pheno$edge.length <- seq_len(nrow(pheno$edge)) / 10

  base_paths <- RERconverge:::allPathsMasterRelativeTT(pheno, master, trees$ap)
  res <- tree2Paths(pheno, trees, binarize = FALSE)
  expect_equal(unname(res), unname(base_paths))
  expect_identical(names(res), colnames(trees$paths))

  bin_tree <- pheno
  bin_tree$edge.length <- ifelse(bin_tree$edge.length > 0, 1, 0)
  bin_res <- tree2Paths(bin_tree, trees, binarize = TRUE)
  expect_true(all(na.omit(bin_res) %in% c(0, 1)))

  cat_tree <- pheno
  cat_tree$edge.length <- seq_len(nrow(cat_tree$edge))
  cat_res <- tree2Paths(cat_tree, trees, categorical = TRUE, binarize = FALSE)
  manual_cat <- {
    prepped <- RERconverge:::prepareTreeForTT(cat_tree, master)
    node_vals <- numeric(length(prepped$tip.label) + prepped$Nnode)
    node_vals[prepped$edge[, 2]] <- prepped$edge.length
    tree_paths <- RERconverge:::allPathsTT(prepped, needIndex = FALSE)
    map <- RERconverge:::matchAllnodesTT(prepped, master)
    tree_paths$nodeId[, 1] <- map[tree_paths$nodeId[, 1], 2]
    tree_paths$nodeId[, 2] <- map[tree_paths$nodeId[, 2], 2]
    ii <- trees$ap$matIndex[cbind(tree_paths$nodeId[, 1], tree_paths$nodeId[, 2])]
    out <- rep(NA_real_, length(trees$ap$dist))
    out[ii] <- node_vals[tree_paths$nodeId[, 1]]
    stats::setNames(out, colnames(trees$paths))
  }
  expect_equal(cat_res, manual_cat)

  trees_no_ap <- trees
  trees_no_ap$ap <- NULL
  expect_equal(tree2Paths(pheno, trees_no_ap), res)
})

test_that("tree2PathsClades and tree2Paths_map match tree2Paths", {
  trees <- build_simple_trees()
  fg_vec <- c("A", "C")
  fg_tree <- foreground2TreeClades(fg_vec, sisters_list = NULL, trees = trees, plotTree = FALSE)
  direct <- tree2Paths(fg_tree, trees)
  expect_equal(tree2PathsClades(fg_tree, trees), direct)

  fake_map <- list(matrix(1))
  expect_equal(tree2Paths_map(fg_tree, fake_map[[1]], trees), direct)
})
