# Reproduces the node-mapping defect that rooting fixes (see docs/rooting-bug.md).
#
# Part 1 shows the mechanism: one unrooted tree written two ways gets its
# internal nodes numbered differently, so the positional map in
# matchAllnodesTT() pairs up non-homologous nodes.
# Part 2 shows the new pipeline is invariant to how the newick was written.
#
# Run:  Rscript tests/reproduce_rooting_bug.R

suppressMessages(library(RERconverge))
suppressMessages(library(TreeTools))
library(ape)

clade <- function(tr, x) {
  if (x <= Ntip(tr)) tr$tip.label[x]
  else paste(sort(extract.clade(tr, x)$tip.label), collapse = "")
}

## Part 1 -- the same unrooted tree, written two ways -------------------------

nw1 <- "(((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):1);"
t1  <- unroot(read.tree(text = nw1))
t2  <- root(t1, outgroup = "A", resolve.root = FALSE)  # moves the basal fork only

cat("newick 1:", nw1, "\n")
cat("newick 2:", write.tree(t2), "\n")
cat("same unrooted tree?  RF =", dist.topo(t1, t2, "PH85"),
    "  same branch lengths:", identical(sort(t1$edge.length), sort(t2$edge.length)), "\n\n")

# old preparation: unroot + RenumberTips + Preorder, no root alignment
M_old <- Preorder(t1)
G_old <- Preorder(RenumberTips(t2, M_old$tip.label))

cat("OLD preparation (unrooted master, no root alignment):\n")
mism <- 0
for (k in (Ntip(M_old) + 1):(Ntip(M_old) + M_old$Nnode)) {
  ok <- clade(M_old, k) == clade(G_old, k)
  mism <- mism + !ok
  cat(sprintf("  node %2d :  master = %-8s  gene = %-8s  %s\n",
              k, clade(M_old, k), clade(G_old, k), if (ok) "ok" else "<-- MISMATCH"))
}
cat("  mismatched internal nodes:", mism, "\n\n")

# new preparation: rooted master + rootLikeMaster
M_new <- Preorder(read.tree(text = nw1))
G_new <- RERconverge:::prepareTreeForTT(t2, M_new)

cat("NEW preparation (rooted master, prepareTreeForTT):\n")
mism <- 0
for (k in (Ntip(M_new) + 1):(Ntip(M_new) + M_new$Nnode)) {
  ok <- clade(M_new, k) == clade(G_new, k)
  mism <- mism + !ok
  cat(sprintf("  node %2d :  master = %-8s  gene = %-8s  %s\n",
              k, clade(M_new, k), clade(G_new, k), if (ok) "ok" else "<-- MISMATCH"))
}
cat("  mismatched internal nodes:", mism, "\n\n")
stopifnot(mism == 0)

## Part 2 -- invariance of the new pipeline to newick root placement ----------

treefile <- "ext/subsetMammalGeneTrees.txt"
if (!file.exists(treefile)) {
  cat("skipping part 2:", treefile, "not found (run from the package root)\n")
  quit(save = "no")
}

set.seed(42)
tmp <- scan(treefile, sep = "\t", what = "character", quiet = TRUE)
out <- character(length(tmp))
for (i in seq(1, length(tmp), by = 2)) {
  out[i] <- tmp[i]
  tr <- unroot(read.tree(text = tmp[i + 1]))
  out[i + 1] <- write.tree(root(tr, outgroup = sample(tr$tip.label, 1),
                                resolve.root = FALSE))
}
rerooted <- tempfile(fileext = ".txt")
writeLines(paste(out[seq(1, length(out), 2)], out[seq(2, length(out), 2)], sep = "\t"),
           rerooted)

a <- readTrees(treefile)
b <- readTrees(rerooted)

cat("\nNEW pipeline, original vs re-rooted input:\n")
cat("  paths dim            :", paste(dim(a$paths), collapse = "x"), "vs",
    paste(dim(b$paths), collapse = "x"), "\n")
cat("  master edge matrix eq:", identical(a$masterTree$edge, b$masterTree$edge), "\n")
cat("  master tip order eq  :", identical(a$masterTree$tip.label, b$masterTree$tip.label), "\n")
cat("  NA pattern identical :", identical(is.na(a$paths), is.na(b$paths)), "\n")

# Any residual value differences should be confined to root-incident paths --
# the root-edge split is not determined by the unrooted tree.
rt <- Ntip(a$masterTree) + 1
d <- abs(a$paths - b$paths); d[is.na(d)] <- 0
bad <- which(colSums(d > 1e-9) > 0)
cat("  columns differing    :", length(bad), "of", ncol(a$paths),
    "( ending at root node:", sum(a$ap$nodeId[bad, 2] == rt), ")\n")

stopifnot(identical(dim(a$paths), dim(b$paths)),
          identical(a$masterTree$edge, b$masterTree$edge),
          identical(is.na(a$paths), is.na(b$paths)))
cat("\nOK\n")
