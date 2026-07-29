simple_tree_lines <- c(
  "geneX\t((A:0.2,B:0.3):0.5,(C:0.4,D:0.6):0.7);",
  "geneY\t((A:0.1,B:0.2):0.4,(C:0.3,D:0.5):0.2);"
)

build_simple_trees <- function(lines = simple_tree_lines) {
  parse_line <- function(line) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    ape::read.tree(text = parts[2])
  }
  trees <- lapply(lines, parse_line)
  master <- TreeTools::Preorder(TreeTools::SortTree(trees[[1]]))
  ap <- RERconverge:::allPathsTT(master)
  treesObj <- list(
    masterTree = master,
    ap = ap,
    matIndex = ap$matIndex,
    maxSp = length(master$tip.label)
  )
  treesObj$paths <- matrix(NA_real_, nrow = length(lines), ncol = length(ap$dist))
  colnames(treesObj$paths) <- RERconverge:::namePathsWSpeciesTT(treesObj)
  class(treesObj) <- c("list", "treesObj")
  treesObj
}

mammal_subset_trees <- function(n = 2) {
  src <- testthat::test_path("..", "..", "ext", "mammal_62_aa_sub.tre")
  if (!file.exists(src)) {
    testthat::skip("mammal_62_aa_sub.tre not available")
  }
  lines <- readLines(src, n)
  parsed <- lapply(lines, function(line) {
    parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
    ape::read.tree(text = parts[2])
  })
  names(parsed) <- vapply(lines, function(line) strsplit(line, "\t", fixed = TRUE)[[1]][1],
                          character(1))
  parsed
}
