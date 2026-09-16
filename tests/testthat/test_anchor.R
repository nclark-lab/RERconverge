# Tests for the anchored path representation, topology status and the fixes that
# came with it. Truth comes from ape (helper-anchor.R), not from the mapping code.

prune_random <- function(t, minKeep) {
  k <- sample(minKeep:ape::Ntip(t), 1)
  ape::unroot(ape::keep.tip(t, sample(t$tip.label, k)))
}

# Gene trees drawn from a master: random prunings, clade removals and a cascade
# that strips the smaller side of the root repeatedly, so that some genes miss a
# whole side of any internal node.
make_genes <- function(M, nGenes, minKeep, regime = "normal") {
  base <- with_lengths(M, regime)
  # each gene has its own rates: an overall scale times per-edge noise, so that
  # residual weighting sees real variance (zero lengths stay zero)
  withRates <- function(t) {
    t$edge.length <- t$edge.length * stats::runif(1, 0.5, 2) * exp(stats::rnorm(nrow(t$edge), 0, 0.3))
    t
  }
  genes <- list(ape::unroot(withRates(base)))
  for (r in seq_len(nGenes)) genes[[length(genes) + 1]] <- prune_random(withRates(base), minKeep)
  mi <- index_master(base)
  internal <- (mi$nT + 2):mi$ntot
  for (v in sample(internal, min(8, length(internal)))) {
    keep <- mi$tips[!mi$desc[v, ]]
    if (length(keep) >= minKeep) genes[[length(genes) + 1]] <- ape::unroot(ape::keep.tip(withRates(base), keep))
  }
  cur <- withRates(base)
  repeat {
    kids <- cur$edge[cur$edge[, 1] == ape::Ntip(cur) + 1, 2]
    sizes <- vapply(kids, function(k) if (k <= ape::Ntip(cur)) 1L else length(ape::extract.clade(cur, k)$tip.label), 0L)
    small <- kids[which.min(sizes)]
    drop <- if (small <= ape::Ntip(cur)) cur$tip.label[small] else ape::extract.clade(cur, small)$tip.label
    if (ape::Ntip(cur) - length(drop) < minKeep) break
    cur <- ape::drop.tip(cur, drop)
    genes[[length(genes) + 1]] <- ape::unroot(cur)
  }
  genes
}

expect_rows_match_truth <- function(tr, newicks, tol = 1e-9) {
  for (i in seq_along(newicks)) {
    tru <- truth_row(tr, newicks[i])
    P <- unname(tr$paths[i, ])  # paths columns carry species names; the truth does not
    expect_identical(!is.na(P), !is.na(tru$v), info = sprintf("row %d filled pattern", i))
    ok <- !is.na(P)
    expect_true(all(abs(P[ok] - tru$v[ok]) <= tol * (1 + abs(tru$v[ok]))),
                info = sprintf("row %d values", i))
  }
}

test_that("anchor master: rooted at a real node, schema from ape", {
  set.seed(11)
  M <- with_lengths(ape::rtree(30))
  genes <- make_genes(M, 60, 12)
  nw <- vapply(genes, to_newick, "")
  # the supplied lengths are kept here, so that rooting can be checked against them
  tr <- read_quiet(write_genes(nw), masterTree = M, masterBranchLengths = "supplied")

  # ape functions misread TreeTools' "preorder" order attribute
  m <- RERconverge:::apeOrder(tr$masterTree)
  expect_equal(as.numeric(ape::dist.topo(ape::unroot(m), ape::unroot(M))), 0)
  sides <- tr$anchor$sides
  expect_setequal(unlist(sides), M$tip.label)
  expect_equal(sum(lengths(sides)), ape::Ntip(M))
  # the root is the anchor node itself: its children are exactly the sides, and
  # no edge was added (no zero-length stem)
  mc <- ape::reorder.phylo(m, "cladewise")
  rootEdges <- which(mc$edge[, 1] == ape::Ntip(mc) + 1)
  childTips <- lapply(mc$edge[rootEdges, 2], function(v) if (v <= ape::Ntip(mc)) mc$tip.label[v] else ape::extract.clade(mc, v)$tip.label)
  expect_equal(length(childTips), length(sides))
  expect_setequal(vapply(childTips, function(x) paste(sort(x), collapse = ","), ""),
                  vapply(sides, function(x) paste(sort(x), collapse = ","), ""))
  expect_equal(nrow(mc$edge), nrow(ape::unroot(M)$edge))
  cm <- ape::cophenetic.phylo(mc); cu <- ape::cophenetic.phylo(M)[rownames(cm), colnames(cm)]
  expect_equal(max(abs(cm - cu)), 0, tolerance = 1e-9)

  # column schema: exactly the ancestor-descendant pairs of the master, one column each
  ep <- expected_pairs(m)
  expect_equal(ncol(tr$paths), nrow(ep))
  cols <- tr$matIndex[ep]
  expect_false(anyNA(cols))
  expect_setequal(cols, seq_len(ncol(tr$paths)))

  # the anchor maximises genes covering every side (brute force with ape)
  u <- ape::unroot(M)
  r <- ape::reorder.phylo(ape::root(u, outgroup = u$tip.label[1], resolve.root = TRUE), "cladewise")
  rmi <- index_master(r)
  pres <- t(vapply(genes, function(g) r$tip.label %in% g$tip.label, logical(ape::Ntip(r))))
  best <- max(vapply((rmi$nT + 2):rmi$ntot, function(v) {
    s <- c(lapply(rmi$kids[[v]], function(k) rmi$desc[k, ]), list(!rmi$desc[v, ]))
    sum(Reduce(`&`, lapply(s, function(side) (pres %*% side) > 0)))
  }, 0))
  expect_equal(tr$anchor$genesAllSides, best)
})

test_that("anchor choice and paths do not depend on how the master is given", {
  set.seed(12)
  M <- with_lengths(ape::rtree(25))
  genes <- make_genes(M, 40, 12)
  f <- write_genes(vapply(genes, to_newick, ""))
  ref <- read_quiet(f, masterTree = M)
  variants <- list(
    unrooted = ape::unroot(M),
    rewritten = ape::read.tree(text = rewrite_newick(M, "tip")),
    edge = ape::read.tree(text = rewrite_newick(M, "edge")))
  for (nm in names(variants)) {
    # an unrooted master warns about midpoint rooting for trait analysis
    tr <- suppressWarnings(read_quiet(f, masterTree = variants[[nm]]))
    expect_identical(ape::write.tree(tr$masterTree), ape::write.tree(ref$masterTree), info = nm)
    expect_equal(tr$paths, ref$paths, tolerance = 1e-12, info = nm)
  }
})

test_that("every filled value equals ape truth: pathological masters and edge lengths", {
  set.seed(13)
  masters <- list(
    rtree40 = ape::rtree(40),
    caterpillar60 = ape::stree(60, "left"),
    balanced32 = ape::stree(32, "balanced"),
    rtree15 = ape::rtree(15))
  nFull <- 0; nDegraded <- 0
  for (mn in names(masters)) {
    for (regime in c("normal", "special", "wide", "equal")) {
      M <- with_lengths(masters[[mn]], regime)
      genes <- make_genes(masters[[mn]], 25, 10, regime)
      nw <- vapply(genes, to_newick, "")
      tr <- read_quiet(write_genes(nw), masterTree = M)
      expect_true(all(tr$treeStatus == "ok"), info = paste(mn, regime))
      full <- vapply(nw, function(x) truth_row(tr, x)$full, NA)
      nFull <- nFull + sum(full); nDegraded <- nDegraded + sum(!full)
      expect_rows_match_truth(tr, nw)
    }
  }
  # the battery must exercise both genes rooted at the anchor node and genes
  # missing a side of it
  expect_gt(nFull, 50)
  expect_gt(nDegraded, 50)
})

test_that("rows do not depend on how each newick is written", {
  set.seed(14)
  M <- with_lengths(ape::rtree(35))
  genes <- make_genes(M, 20, 10)
  pick <- sample(seq_along(genes), 10)
  nw <- vapply(genes, to_newick, "")
  extra <- unlist(lapply(pick, function(j) vapply(c("tip", "edge", "unrooted", "ladder"),
                                                  function(h) rewrite_newick(genes[[j]], h), "")))
  tr <- read_quiet(write_genes(c(nw, extra)), masterTree = M)
  expect_true(all(tr$treeStatus == "ok"))
  k <- length(nw)
  for (a in seq_along(pick)) {
    for (h in 1:4) {
      k <- k + 1
      expect_identical(is.na(tr$paths[k, ]), is.na(tr$paths[pick[a], ]))
      expect_equal(tr$paths[k, ], tr$paths[pick[a], ], tolerance = 1e-12)
    }
  }
})

test_that("pruned-pair branch vectors read from paths equal direct pruning", {
  set.seed(15)
  M <- with_lengths(ape::rtree(45), "special")
  genes <- make_genes(M, 50, 12, "special")
  nw <- vapply(genes, to_newick, "")
  tr <- read_quiet(write_genes(nw), masterTree = M)
  for (p in 1:60) {
    ij <- sample(seq_along(nw), 2)
    C <- intersect(genes[[ij[1]]]$tip.label, genes[[ij[2]]]$tip.label)
    if (length(C) < 4) next
    for (k in ij) {
      want <- pruned_edges(nw[k], C)
      got <- row_edges(tr, k, C)
      expect_setequal(names(got), names(want))
      expect_equal(unname(got[names(want)]), unname(want), tolerance = 1e-9)
    }
  }
})

test_that("real-node values do not move when species are dropped", {
  set.seed(16)
  M <- with_lengths(ape::rtree(40))
  base <- ape::unroot(with_lengths(M))
  subsets <- replicate(30, sample(base$tip.label, sample(12:39, 1)), simplify = FALSE)
  nw <- c(to_newick(base), vapply(subsets, function(s) to_newick(ape::unroot(ape::keep.tip(base, s))), ""))
  tr <- read_quiet(write_genes(nw), masterTree = M)
  b <- truth_row(tr, nw[1])
  expect_true(b$full)
  for (i in 2:length(nw)) {
    t_i <- truth_row(tr, nw[i])
    both <- !is.na(tr$paths[1, ]) & !is.na(tr$paths[i, ])
    # a column may move only if it ends at this gene's root and that root is a
    # point on an edge (the gene misses a side of the anchor)
    movable <- if (t_i$full) rep(FALSE, length(both)) else t_i$atRoot
    stable <- both & !movable
    expect_equal(tr$paths[i, stable], tr$paths[1, stable], tolerance = 1e-9, info = sprintf("subset %d", i))
  }
})

test_that("discordance is decided by topology and decoupled from path lookup", {
  set.seed(17)
  M <- with_lengths(ape::rtree(30))
  base <- ape::unroot(with_lengths(M))
  good <- make_genes(M, 15, 12)
  # discordant: swap two tips from opposite sides of the root
  mc <- ape::reorder.phylo(M, "cladewise")
  rk <- mc$edge[mc$edge[, 1] == ape::Ntip(mc) + 1, 2]
  sideOf <- function(v) if (v <= ape::Ntip(mc)) mc$tip.label[v] else ape::extract.clade(mc, v)$tip.label
  swap <- base
  x <- sideOf(rk[1])[1]; y <- sideOf(rk[2])[1]
  swap$tip.label[match(c(x, y), swap$tip.label)] <- c(y, x)
  # discordant: nearest-neighbour interchange on an internal edge
  nni <- ape::unroot(phangorn::nni(base)[[1]])
  # unresolved: collapse one internal edge
  poly <- base
  inner <- which(poly$edge[, 2] > ape::Ntip(poly))
  poly$edge.length[inner[1]] <- 0
  poly <- ape::di2multi(poly)
  # concordant genes written with every layout, including degraded anchors
  layouts <- unlist(lapply(good[1:5], function(g) vapply(c("tip", "edge", "unrooted", "ladder"), function(h) rewrite_newick(g, h), "")))

  nw <- c(vapply(good, to_newick, ""), to_newick(swap), to_newick(nni), to_newick(poly), layouts)
  geneNames <- sprintf("g%03d", seq_along(nw))
  flagged <- length(good) + 1:3
  expect_no_error(tr <- read_quiet(write_genes(nw), masterTree = M))

  # trees that do not match the master are dropped, with their reason recorded
  expect_equal(tr$numTrees, length(nw) - length(flagged))
  expect_setequal(tr$dropped$gene, geneNames[flagged])
  expect_identical(tr$dropped$reason[match(geneNames[flagged], tr$dropped$gene)],
                   c("discordant", "discordant", "unresolved"))
  expect_false(any(geneNames[flagged] %in% names(tr$trees)))

  # everything kept has paths, and they match the ape truth
  expect_true(all(tr$treeStatus == "ok"))
  expect_true(all(rowSums(!is.na(tr$paths)) > 0))
  expect_rows_match_truth(tr, nw[-flagged][seq_len(length(good))])
})

test_that("a failed path lookup on a concordant tree is an error, not a discordance flag", {
  set.seed(18)
  M <- with_lengths(ape::rtree(20))
  genes <- make_genes(M, 10, 12)
  # few genes here: keep the supplied lengths instead of estimating
  tr <- read_quiet(write_genes(vapply(genes, to_newick, "")), masterTree = M,
                   masterBranchLengths = "supplied")
  g <- genes[[1]]
  # concordant, but rooted far from the master's root and never prepared
  bad <- TreeTools::Preorder(TreeTools::RenumberTips(
    ape::root(g, outgroup = g$tip.label[3], resolve.root = TRUE), tr$masterTree$tip.label))
  expect_identical(RERconverge:::treeTopologyStatus(bad, tr$masterTree), "ok")
  expect_error(RERconverge:::allPathsMasterRelativeTT(bad, tr$masterTree, tr$ap, 1, check_concordance = FALSE),
               "internal error")
})

test_that("small trees are dropped at input", {
  set.seed(19)
  M <- with_lengths(ape::rtree(30))
  genes <- make_genes(M, 20, 12)
  tiny <- ape::unroot(ape::keep.tip(with_lengths(M), M$tip.label[1:6]))
  nw <- c(vapply(genes, to_newick, ""), to_newick(tiny))
  expect_message(tr <- readTrees(write_genes(nw), masterTree = M), "Dropped 1 tree")
  expect_equal(tr$numTrees, length(genes))
  expect_equal(nrow(tr$paths), length(genes))
  # the dropped tree is recorded, with its reason, and is not among the trees
  expect_equal(nrow(tr$dropped), 1)
  expect_identical(tr$dropped$reason, "few_species")
  expect_false(tr$dropped$gene %in% names(tr$trees))
  tr5 <- read_quiet(write_genes(nw), masterTree = M, minTreeSpecies = 5)
  expect_equal(tr5$numTrees, length(genes) + 1)
  expect_equal(nrow(tr5$dropped), 0)
})

test_that("master species absent from every gene do not break naming or residuals", {
  set.seed(20)
  M <- with_lengths(ape::rtree(40))
  sub <- ape::keep.tip(M, M$tip.label[1:34])
  genes <- make_genes(sub, 150, 22)
  tr <- read_quiet(write_genes(vapply(genes, to_newick, "")), masterTree = M)
  named <- unique(colnames(tr$paths)[colnames(tr$paths) != ""])
  expect_setequal(named, M$tip.label)
  expect_no_error(rer <- suppressMessages(suppressWarnings(getAllResiduals(tr))))
  expect_equal(dim(rer), dim(tr$paths))
})

test_that("residual regression uses each gene's own branches, with and without useSpecies", {
  set.seed(21)
  M <- with_lengths(ape::rtree(40))
  genes <- make_genes(M, 180, 24)
  nw <- vapply(genes, to_newick, "")
  tr <- read_quiet(write_genes(nw), masterTree = M)

  # imputation warns about sparse rows; that is expected for pruned genes
  tr2 <- suppressMessages(suppressWarnings(transformPaths(tr)))
  res <- suppressMessages(suppressWarnings(coreGetResiduals(tr2)))
  checked <- 0
  for (i in seq_along(nw)) {
    idx <- res$index[[i]]
    if (!length(idx)) next
    own <- truth_row(tr, nw[i])$edgeCols
    expect_true(all(idx %in% own), info = sprintf("gene %d", i))
    checked <- checked + 1
  }
  expect_gt(checked, 100)

  useSp <- sample(M$tip.label, 30)
  resU <- suppressMessages(suppressWarnings(coreGetResiduals(tr2, useSpecies = useSp)))
  checked <- 0
  for (i in seq_along(nw)) {
    idx <- resU$index[[i]]
    if (!length(idx)) next
    keep <- intersect(genes[[i]]$tip.label, useSp)
    own <- truth_row(tr, to_newick(ape::unroot(ape::keep.tip(genes[[i]], keep))))$edgeCols
    expect_true(all(idx %in% own), info = sprintf("gene %d with useSpecies", i))
    checked <- checked + 1
  }
  expect_gt(checked, 20)
})

test_that("master branch lengths come from the data unless the supplied tree is used explicitly", {
  set.seed(22)
  M <- with_lengths(ape::rtree(30))
  genes <- make_genes(M, 60, 15)
  f <- write_genes(vapply(genes, to_newick, ""))

  # default: estimated from the genes, not taken from the supplied tree
  tr <- read_quiet(f, masterTree = M)
  expect_true(all(is.finite(tr$masterTree$edge.length)))
  expect_true(all(tr$masterTree$edge.length >= 0))
  supplied <- ape::cophenetic.phylo(M)
  estimated <- ape::cophenetic.phylo(RERconverge:::apeOrder(tr$masterTree))
  expect_false(isTRUE(all.equal(estimated, supplied[rownames(estimated), colnames(estimated)])))

  # explicit override keeps the supplied lengths exactly
  trS <- read_quiet(f, masterTree = M, masterBranchLengths = "supplied")
  kept <- ape::cophenetic.phylo(RERconverge:::apeOrder(trS$masterTree))
  expect_equal(kept, supplied[rownames(kept), colnames(kept)], tolerance = 1e-9)

  # the override needs a tree that has branch lengths
  noLen <- M; noLen$edge.length <- NULL
  expect_error(read_quiet(f, masterTree = noLen, masterBranchLengths = "supplied"),
               "branch lengths")

  # minSpecs restricts which genes contribute
  trM <- read_quiet(f, masterTree = M, minSpecs = 25, minTreesAll = 5)
  expect_true(all(is.finite(trM$masterTree$edge.length)))
  expect_false(isTRUE(all.equal(trM$masterTree$edge.length, tr$masterTree$edge.length)))

  # too few usable genes is an error, not a silently unusable master
  expect_error(read_quiet(f, masterTree = M, minTreesAll = 1000), "cannot estimate")

  # the deprecated argument still selects the supplied lengths
  expect_warning(trD <- readTrees(f, masterTree = M, reestimateBranches = FALSE), "deprecated")
  kept2 <- ape::cophenetic.phylo(RERconverge:::apeOrder(trD$masterTree))
  expect_equal(kept2, supplied[rownames(kept2), colnames(kept2)], tolerance = 1e-9)
})

test_that("anchor = 'root' keeps the supplied root and requires one", {
  set.seed(23)
  M <- with_lengths(ape::rtree(25))
  genes <- make_genes(M, 30, 12)
  nw <- vapply(genes, to_newick, "")
  f <- write_genes(nw)
  tr <- read_quiet(f, masterTree = M, anchor = "root")
  expect_null(tr$anchor)
  mc <- ape::reorder.phylo(tr$masterTree, "cladewise")
  rk <- mc$edge[mc$edge[, 1] == ape::Ntip(mc) + 1, 2]
  side <- if (rk[1] <= ape::Ntip(mc)) mc$tip.label[rk[1]] else ape::extract.clade(mc, rk[1])$tip.label
  expect_true(RERconverge:::isRootedOn(M, side))
  expect_rows_match_truth(tr, nw)
  expect_error(read_quiet(f, masterTree = ape::unroot(M), anchor = "root"), "rooted")
})

test_that("rootLikeMaster roots four-species layouts that TreeTools::RootTree leaves unrooted", {
  master <- TreeTools::Preorder(ape::read.tree(text = "((a:1,b:1):1,(c:1,d:1):1);"))
  for (nw in c("(a:1,b:2,(c:3,d:4):5);", "(c:3,d:4,(a:1,b:2):5);", "((a:1,b:2):5,c:3,d:4);")) {
    g <- RERconverge:::rootLikeMaster(ape::read.tree(text = nw), master)
    expect_true(RERconverge:::isRootedOn(g, c("a", "b")), info = nw)
  }
})

test_that("anchor choice does not depend on newick child order (ties)", {
  # four clades of three species; every internal node covers all genes equally,
  # so the anchor is decided by the tie-break alone
  expand <- function(x) sprintf("(%s1:1,%s2:1,%s3:1):1", x, x, x)
  nwA <- sprintf("((%s,%s):1,(%s,%s):1);", expand("a"), expand("z"), expand("b"), expand("y"))
  nwB <- sprintf("((%s,%s):1,(%s,%s):1);", expand("z"), expand("a"), expand("y"), expand("b"))
  MA <- ape::read.tree(text = nwA); MB <- ape::read.tree(text = nwB)
  genes <- rep(list(ape::unroot(MA)), 3)
  f <- write_genes(vapply(genes, to_newick, ""))
  # only a handful of genes here: keep the supplied lengths instead of estimating
  trA <- read_quiet(f, masterTree = MA, masterBranchLengths = "supplied")
  trB <- read_quiet(f, masterTree = MB, masterBranchLengths = "supplied")
  key <- function(s) sort(vapply(s, function(x) paste(sort(x), collapse = ","), ""))
  expect_identical(key(trA$anchor$sides), key(trB$anchor$sides))
  expect_identical(ape::write.tree(trA$masterTree), ape::write.tree(trB$masterTree))
  expect_equal(trA$paths, trB$paths)
})

test_that("multifurcating master: pruning a side of the root maps correctly", {
  set.seed(24)
  # a four-way root and an internal polytomy; genes keep the polytomies
  M <- ape::read.tree(text = paste0(
    "((a1:1,a2:1,a3:1):2,((b1:1,b2:1):1,(b3:1,b4:1):1):3,(c1:1,c2:1,(c3:1,c4:1):1):4,",
    "(d1:1,(d2:1,d3:1):1):5);"))
  genes <- make_genes(M, 30, 10)
  # add genes without each whole root side
  mc <- ape::reorder.phylo(M, "cladewise")
  for (k in mc$edge[mc$edge[, 1] == ape::Ntip(mc) + 1, 2]) {
    drop <- if (k <= ape::Ntip(mc)) mc$tip.label[k] else ape::extract.clade(mc, k)$tip.label
    if (ape::Ntip(M) - length(drop) >= 10) genes[[length(genes) + 1]] <- ape::unroot(ape::drop.tip(M, drop))
  }
  nw <- vapply(genes, to_newick, "")
  # this master is unrooted, so trait analysis warns about midpoint rooting
  expect_no_error(tr <- suppressWarnings(read_quiet(write_genes(nw), masterTree = M)))
  expect_true(all(tr$treeStatus == "ok"))
  expect_rows_match_truth(tr, nw)
})

test_that("kept master vertices are correct when a side of a multifurcating root is dropped", {
  # TreeTools::KeptVerts/KeepTip suppress the first child instead of the root here
  m <- TreeTools::Preorder(ape::read.tree(text = "((a:1,b:2):3,((c:1,d:2):3,g:4):5,(e:1,f:2):6);"))
  mc <- RERconverge:::apeOrder(m)
  clade <- function(t, v) if (v <= ape::Ntip(t)) t$tip.label[v] else sort(ape::extract.clade(t, v)$tip.label)
  for (keep in list(c("a", "b", "c", "d", "g"), c("a", "b", "e", "f"), c("c", "d", "g", "e"), c("c", "d", "g"))) {
    tree <- ape::keep.tip(mc, keep)
    map <- RERconverge:::matchAllnodesTT(tree, m)
    kept <- map[, 2]
    expected <- which(vapply(seq_len(ape::Ntip(mc) + mc$Nnode), function(v) {
      if (v <= ape::Ntip(mc)) return(mc$tip.label[v] %in% keep)
      kids <- mc$edge[mc$edge[, 1] == v, 2]
      sum(vapply(kids, function(k) any(clade(mc, k) %in% keep), NA)) >= 2
    }, NA))
    expect_setequal(kept, expected)
  }
})

test_that("RER tree export handles genes without paths", {
  set.seed(25)
  M <- with_lengths(ape::rtree(20))
  genes <- make_genes(M, 12, 12)
  base <- ape::unroot(with_lengths(M))
  # discordant: swap two tips from opposite sides of the root
  mc <- ape::reorder.phylo(M, "cladewise")
  rk <- mc$edge[mc$edge[, 1] == ape::Ntip(mc) + 1, 2]
  sideOf <- function(v) if (v <= ape::Ntip(mc)) mc$tip.label[v] else ape::extract.clade(mc, v)$tip.label
  swap <- base
  x <- sideOf(rk[1])[1]; y <- sideOf(rk[2])[1]
  swap$tip.label[match(c(x, y), swap$tip.label)] <- c(y, x)
  # a species the master does not have
  other <- base
  other$tip.label[1] <- "OTHER"
  nw <- c(vapply(genes, to_newick, ""), to_newick(swap), to_newick(other))
  tr0 <- read_quiet(write_genes(nw), masterTree = M, masterBranchLengths = "supplied")
  # readTrees drops those two trees
  expect_equal(tr0$numTrees, length(genes))
  expect_setequal(tr0$dropped$reason, c("discordant", "species_not_in_master"))

  # Export must still cope with a treesObj that carries trees without paths:
  # objects saved before such trees were dropped, or built by other code.
  tr <- tr0
  flagged <- c(2L, 5L)
  for (i in flagged) {
    tr$treeStatus[i] <- "discordant"
    tr$paths[i, ] <- NA
    tr$trees[[i]] <- ape::unroot(ape::read.tree(text = nw[i]))  # unprepared, as readTrees left it
  }

  r <- tr$paths
  rownames(r) <- names(tr$trees)

  for (i in flagged) {
    expect_warning(t_i <- returnRersAsTree(tr, r, i, plot = FALSE), "has no RERs")
    expect_s3_class(t_i, "phylo")
    expect_true(all(is.na(t_i$edge.length)))
    expect_setequal(t_i$tip.label, tr$trees[[i]]$tip.label)
    # by name as well as by index
    expect_warning(returnRersAsTree(tr, r, names(tr$trees)[i], plot = FALSE), "has no RERs")
  }

  # mapped genes keep their own branch values
  for (i in setdiff(seq_along(genes), flagged)[1:5]) {
    t_i <- expect_silent(returnRersAsTree(tr, r, i, plot = FALSE))
    own <- r[i, truth_row(tr, nw[i])$edgeCols]
    expect_equal(sort(unname(t_i$edge.length)), sort(unname(own)))
  }

  expect_warning(all <- returnRersAsTreesAll(tr, r), "2 genes have no RERs")
  expect_s3_class(all, "multiPhylo")
  expect_equal(length(all), tr$numTrees)
  expect_true(all(vapply(all[flagged], function(t) all(is.na(t$edge.length)), NA)))
  expect_false(any(vapply(all[-flagged], function(t) all(is.na(t$edge.length)), NA)))

  expect_warning(nwk <- returnRersAsNewickStrings(tr, r), "2 genes have no RERs")
  expect_equal(length(nwk), tr$numTrees)
  expect_false(anyNA(nwk))
})

test_that("the master's tip numbers are contiguous within every clade, and violations are caught", {
  # matchAllnodesTT pairs gene nodes with surviving master nodes by position. That
  # is valid only if pruning never reorders sibling subtrees under Preorder's
  # lowest-numbered-leaf rule, i.e. if every clade's tips are numbered contiguously.
  set.seed(26)
  for (M in list(with_lengths(ape::rtree(40)), with_lengths(ape::stree(30, "left")),
                 with_lengths(ape::stree(32, "balanced")))) {
    genes <- make_genes(M, 12, 12)
    tr <- read_quiet(write_genes(vapply(genes, to_newick, "")), masterTree = M)
    expect_silent(RERconverge:::assertContiguousTips(tr$masterTree))
    tr2 <- read_quiet(write_genes(vapply(genes, to_newick, "")), masterTree = M, anchor = "root")
    expect_silent(RERconverge:::assertContiguousTips(tr2$masterTree))
  }
  # clade {a, c} holds tips 1 and 3
  bad <- structure(list(edge = rbind(c(5L, 6L), c(6L, 1L), c(6L, 3L), c(5L, 7L), c(7L, 2L), c(7L, 4L)),
                        tip.label = c("a", "b", "c", "d"), Nnode = 3L), class = "phylo")
  expect_error(RERconverge:::assertContiguousTips(bad), "not numbered contiguously")
})

test_that("the master keeps a biological rooting for trait analysis", {
  set.seed(27)
  M <- with_lengths(ape::rtree(30))
  genes <- make_genes(M, 40, 12)
  f <- write_genes(vapply(genes, to_newick, ""))
  tr <- read_quiet(f, masterTree = M)

  R <- RERconverge:::apeOrder(tr$masterTreeRooted)
  expect_true(ape::is.rooted(R))
  # same tree as the anchored master, rooted where the supplied tree was rooted
  expect_equal(as.numeric(ape::dist.topo(ape::unroot(R), ape::unroot(M))), 0)
  side <- RERconverge:::rootSides(RERconverge:::apeOrder(M))[[1]]
  expect_true(RERconverge:::isRootedOn(R, side))
  # with the master's estimated branch lengths, not the supplied tree's
  expect_equal(sum(R$edge.length), sum(tr$masterTree$edge.length), tolerance = 1e-9)
  cm <- ape::cophenetic.phylo(R)
  ca <- ape::cophenetic.phylo(RERconverge:::apeOrder(tr$masterTree))[rownames(cm), colnames(cm)]
  expect_equal(max(abs(cm - ca)), 0, tolerance = 1e-9)

  # no biological root available: midpoint, with a warning
  expect_warning(trU <- suppressMessages(readTrees(f, masterTree = ape::unroot(M))), "midpoint")
  expect_true(ape::is.rooted(RERconverge:::apeOrder(trU$masterTreeRooted)))
  expect_equal(as.numeric(ape::dist.topo(ape::unroot(RERconverge:::apeOrder(trU$masterTreeRooted)),
                                         ape::unroot(M))), 0)
})

test_that("char2Paths orients trait change by the biological root, not the anchor", {
  set.seed(28)
  M <- with_lengths(ape::rtree(24))          # its rooting is the biological one
  genes <- make_genes(M, 60, 12)
  tr <- read_quiet(write_genes(vapply(genes, to_newick, "")), masterTree = M)
  tips <- stats::setNames(stats::rnorm(ape::Ntip(M)), M$tip.label)
  p <- char2Paths(tips, tr)
  expect_equal(length(p), ncol(tr$paths))
  expect_true(any(!is.na(p)))

  # independent ape computation of the states on the biologically rooted master
  R <- RERconverge:::apeOrder(tr$masterTreeRooted)
  A <- RERconverge:::apeOrder(tr$masterTree)
  fa <- phytools::fastAnc(R, tips[R$tip.label])
  state <- c(tips[R$tip.label], fa)
  cladeOf <- function(t, v) if (v <= ape::Ntip(t)) t$tip.label[v] else ape::extract.clade(t, v)$tip.label
  key <- function(x) { x <- sort(x); if (sort(A$tip.label)[1] %in% x) x <- setdiff(A$tip.label, x); paste(sort(x), collapse = ",") }
  # anchored single-edge columns, by bipartition
  colOf <- stats::setNames(tr$matIndex[cbind(A$edge[, 2], A$edge[, 1])],
                           vapply(A$edge[, 2], function(v) key(cladeOf(A, v)), ""))
  rootKids <- R$edge[R$edge[, 1] == ape::Ntip(R) + 1L, 2]

  checked <- 0
  for (k in seq_len(nrow(R$edge))) {
    pR <- R$edge[k, 1]; cR <- R$edge[k, 2]
    if (pR == ape::Ntip(R) + 1L) next               # the two halves of the root branch
    col <- colOf[[key(cladeOf(R, cR))]]
    # biological direction: descendant minus ancestor, whatever the anchored one is
    expect_equal(unname(p[col]), unname(state[cR] - state[pR]), tolerance = 1e-8)
    checked <- checked + 1
  }
  expect_gt(checked, 30)

  # the branch containing the biological root: difference across the whole branch
  colRoot <- colOf[[key(cladeOf(R, rootKids[1]))]]
  expect_equal(abs(unname(p[colRoot])),
               abs(unname(state[rootKids[1]] - state[rootKids[2]])), tolerance = 1e-8)
})

test_that("concordant_trees accepts TreeTools-preordered trees", {
  m <- TreeTools::Preorder(ape::read.tree(text = "((a:1,b:1):1,((c:1,d:1):1,(e:1,f:1):1):1);"))
  g <- TreeTools::Preorder(ape::read.tree(text = "((f:1,e:1):1,(d:1,c:1):1,(b:1,a:1):1);"))
  expect_true(concordant_trees(g, m))
})
