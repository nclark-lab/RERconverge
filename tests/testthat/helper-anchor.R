# Shared helpers for the anchor tests. Truth is computed with ape only, never with
# the TreeTools-based mapping under test.

write_genes <- function(newicks, names = sprintf("g%03d", seq_along(newicks))) {
  f <- tempfile(fileext = ".txt")
  writeLines(paste0(names, "\t", newicks), f)
  f
}

read_quiet <- function(f, ...) {
  suppressMessages(readTrees(f, ...))
}

# Full-precision newick, optionally with shuffled child order, built bottom-up
# (recursion overflows the C stack on deep caterpillars).
to_newick <- function(t, shuffle = FALSE) {
  n <- ape::Ntip(t); N <- n + t$Nnode
  lab <- character(N); lab[1:n] <- t$tip.label
  len <- character(N); len[t$edge[, 2]] <- sprintf("%.17g", t$edge.length)
  kids <- split(t$edge[, 2], factor(t$edge[, 1], levels = 1:N))
  depth <- integer(N)
  e <- ape::reorder.phylo(t, "cladewise")$edge
  for (k in seq_len(nrow(e))) depth[e[k, 2]] <- depth[e[k, 1]] + 1L
  internal <- (n + 1):N
  for (v in internal[order(-depth[internal])]) {
    ks <- kids[[v]]
    if (shuffle) ks <- ks[sample.int(length(ks))]
    lab[v] <- paste0("(", paste0(lab[ks], ":", len[ks], collapse = ","), ")")
  }
  paste0(lab[n + 1], ";")
}

# The same unrooted tree written four ways.
rewrite_newick <- function(t, how) {
  u <- ape::unroot(t)
  x <- switch(how,
    tip = ape::root(u, outgroup = sample(u$tip.label, 1), resolve.root = TRUE),
    edge = {
      r <- ape::root(u, outgroup = u$tip.label[1], resolve.root = TRUE)
      inner <- (ape::Ntip(r) + 2):(ape::Ntip(r) + r$Nnode)
      sizes <- vapply(inner, function(v) length(ape::extract.clade(r, v)$tip.label), 0L)
      ok <- inner[sizes >= 2 & sizes <= ape::Ntip(r) - 2]
      if (length(ok)) ape::root(r, outgroup = ape::extract.clade(r, ok[sample.int(length(ok), 1)])$tip.label,
                                resolve.root = TRUE) else r
    },
    unrooted = u,
    ladder = ape::ladderize(ape::root(u, outgroup = sample(u$tip.label, 1), resolve.root = TRUE)))
  to_newick(x, shuffle = how != "ladder")
}

with_lengths <- function(t, regime = c("normal", "special", "wide", "equal")) {
  regime <- match.arg(regime)
  ne <- nrow(t$edge)
  t$edge.length <- switch(regime,
    normal = stats::rexp(ne, 10),
    special = { x <- stats::rexp(ne, 10); k <- sample(ne, ceiling(0.15 * ne)); x[k] <- rep(c(0, 1e-8), length.out = length(k)); x },
    wide = 10^stats::runif(ne, -12, 1),
    equal = rep(1, ne))
  t
}

# Structure of a rooted master: tips, children, descendant tips, MRCA table.
index_master <- function(m) {
  m <- ape::reorder.phylo(m, "cladewise")
  nT <- ape::Ntip(m); ntot <- nT + m$Nnode; mroot <- nT + 1L
  kids <- split(m$edge[, 2], factor(m$edge[, 1], levels = 1:ntot))
  parent <- integer(ntot); parent[m$edge[, 2]] <- m$edge[, 1]
  anc <- lapply(1:ntot, function(v) { a <- v; while (v != mroot) { v <- parent[v]; a <- c(a, v) }; a })
  desc <- matrix(FALSE, ntot, nT)
  for (t in 1:nT) desc[anc[[t]], t] <- TRUE
  ancMat <- matrix(FALSE, ntot, ntot)
  for (v in 1:ntot) ancMat[v, anc[[v]]] <- TRUE
  mrca <- matrix(0L, ntot, ntot)
  for (a in 1:ntot) mrca[a, ] <- anc[[a]][max.col(ancMat[, anc[[a]], drop = FALSE], "first")]
  list(tips = m$tip.label, nT = nT, ntot = ntot, kids = kids, anc = anc,
       depth = lengths(anc) - 1L, desc = desc, mrca = mrca)
}

# Every ancestor-descendant pair of the master, from ape alone.
expected_pairs <- function(m) {
  mi <- index_master(m)
  pairs <- do.call(rbind, lapply(seq_len(mi$ntot), function(d) {
    a <- mi$anc[[d]][-1]
    if (length(a)) cbind(d, a) else NULL
  }))
  pairs[order(pairs[, 1], pairs[, 2]), , drop = FALSE]
}

column_pairs <- function(tr) {
  pos <- which(tr$matIndex > 0, arr.ind = TRUE)
  cp <- matrix(NA_integer_, ncol(tr$paths), 2)
  cp[tr$matIndex[pos], ] <- pos
  cp
}

# The row readTrees should produce for a concordant gene, computed with ape:
# root the gene at the master's root restricted to its species, split the root
# edge as documented (whole edge on anchor side A when the gene covers every side
# of the anchor, evenly otherwise), map nodes by clade, and read node distances.
truth_row <- function(tr, newick) {
  m <- tr$masterTree
  mi <- index_master(m)
  g <- ape::unroot(ape::read.tree(text = newick))
  S <- match(g$tip.label, mi$tips)
  cnt <- rowSums(mi$desc[, S, drop = FALSE])
  proot <- which(cnt == length(S)); proot <- proot[which.max(mi$depth[proot])]
  ks <- mi$kids[[proot]]; ks <- ks[cnt[ks] > 0]
  sideOf <- function(k) mi$tips[mi$desc[k, ] & mi$tips %in% g$tip.label]
  # The gene is rooted where the master root falls after pruning: at a node when
  # it keeps three or more sides (every value is then a real distance), otherwise
  # on the edge it has become, split evenly.
  full <- length(ks) >= 3
  if (full) {
    small <- sideOf(ks[which.min(vapply(ks, function(k) length(sideOf(k)), 0L))])
    g1 <- ape::reorder.phylo(ape::root(g, outgroup = setdiff(g$tip.label, small)[1], resolve.root = TRUE), "cladewise")
    anc <- if (length(small) == 1) match(small, g1$tip.label) else ape::getMRCA(g1, small)
    g2 <- ape::root(g1, node = g1$edge[g1$edge[, 2] == anc, 1], resolve.root = FALSE)
    stopifnot(sum(g2$edge[, 1] == ape::Ntip(g2) + 1L) == length(ks))
  } else {
    a <- sideOf(ks[1])
    b <- setdiff(g$tip.label, a)
    og <- if (length(a) <= length(b)) a else b
    g1 <- ape::root(g, outgroup = setdiff(g$tip.label, og)[1], resolve.root = TRUE)
    g2 <- ape::root(g1, outgroup = og, resolve.root = TRUE)
    stopifnot(ape::is.monophyletic(g2, og))
    rootEdges <- which(g2$edge[, 1] == ape::Ntip(g2) + 1L)
    g2$edge.length[rootEdges] <- sum(g2$edge.length[rootEdges]) / 2
  }
  n <- ape::Ntip(g2)
  gN <- n + g2$Nnode
  gm <- integer(gN); gm[1:n] <- match(g2$tip.label, mi$tips)
  e <- ape::reorder.phylo(g2, "postorder")$edge
  for (k in seq_len(nrow(e))) {
    p <- e[k, 1]; c <- e[k, 2]
    gm[p] <- if (gm[p] == 0L) gm[c] else mi$mrca[gm[p], gm[c]]
  }
  inv <- integer(mi$ntot); inv[gm] <- seq_len(gN)
  nk <- vapply(mi$kids, function(k) sum(cnt[k] > 0), 0L)
  kept <- c(S, which(nk >= 2 & seq_len(mi$ntot) > mi$nT))
  cp <- column_pairs(tr)
  E <- cp[, 1] %in% kept & cp[, 2] %in% kept
  stopifnot(all(inv[kept] > 0))
  D <- ape::dist.nodes(g2)
  v <- rep(NA_real_, nrow(cp))
  v[E] <- D[cbind(inv[cp[E, 1]], inv[cp[E, 2]])]
  atRoot <- cp[, 1] == proot | cp[, 2] == proot
  # the gene's own edges, as (descendant, ancestor) columns
  edgeCols <- tr$matIndex[cbind(gm[g2$edge[, 2]], gm[g2$edge[, 1]])]
  list(v = v, full = full, proot = proot, atRoot = atRoot, rootKids = ks, cp = cp,
       edgeCols = edgeCols)
}

# Branch lengths of a tree pruned to species C, keyed by bipartition (ape only).
pruned_edges <- function(newick, C) {
  C <- sort(C)
  g <- ape::unroot(ape::keep.tip(ape::unroot(ape::read.tree(text = newick)), C))
  n <- ape::Ntip(g)
  below <- vector("list", n + g$Nnode)
  for (k in 1:n) below[[k]] <- g$tip.label[k]
  po <- ape::reorder.phylo(g, "postorder")
  for (k in seq_len(nrow(po$edge))) below[[po$edge[k, 1]]] <- c(below[[po$edge[k, 1]]], below[[po$edge[k, 2]]])
  key <- function(x) { x <- sort(x); if (C[1] %in% x) x <- setdiff(C, x); paste(x, collapse = ",") }
  stats::setNames(po$edge.length, vapply(po$edge[, 2], function(v) key(below[[v]]), ""))
}

# The same edges read from a paths row: one column per edge of the master reduced
# to C, plus the sum of the two columns at the reduced root when it has two
# children.
row_edges <- function(tr, i, C) {
  C <- sort(C)
  mi <- index_master(tr$masterTree)
  S <- match(C, mi$tips)
  cnt <- rowSums(mi$desc[, S, drop = FALSE])
  nk <- vapply(mi$kids, function(k) sum(cnt[k] > 0), 0L)
  kept <- c(S, which(nk >= 2 & seq_len(mi$ntot) > mi$nT))
  proot <- kept[which.max(ifelse(cnt[kept] == length(S), mi$depth[kept], -1))]
  parentKept <- function(v) { a <- mi$anc[[v]][-1]; a[a %in% kept][1] }
  nonroot <- setdiff(kept, proot)
  pk <- vapply(nonroot, parentKept, 0L)
  val <- tr$paths[i, tr$matIndex[cbind(nonroot, pk)]]
  key <- function(v) { x <- sort(mi$tips[mi$desc[v, ] & mi$tips %in% C]); if (C[1] %in% x) x <- setdiff(C, x); paste(x, collapse = ",") }
  keys <- vapply(nonroot, key, "")
  out <- stats::setNames(val, keys)
  atRoot <- pk == proot
  if (sum(atRoot) == 2) {
    out <- c(out[!atRoot], stats::setNames(sum(val[atRoot]), keys[atRoot][1]))
  }
  out
}
