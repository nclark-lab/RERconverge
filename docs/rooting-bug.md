# Why all trees are now rooted

Background note for the `feature-rooted-trees` branch (PR #127). Describes the
defect in the path-indexing machinery that rooting fixes, how to reproduce it,
and what it does and does not fix.

Reproducer: `tests/reproduce_rooting_bug.R`.

**Scope.** This rarely bites in practice. Newick strings are normally produced
by one pipeline in one consistent way, so all the trees in a file carry their
basal fork in the same place, the positional map happens to be correct, and
nothing goes wrong — the bundled mammal dataset is exactly such a case (section
4). The defect surfaces when that convention is not uniform: trees assembled
from more than one source or tool, files edited or re-exported through software
that re-roots on write, phenotype/trait trees built separately from the gene
trees, or a master tree supplied from elsewhere. Because the failure is silent
when it happens (section 3c), correctness should not rest on the input
convention holding.

---

## 1. The indexing scheme

RERconverge represents every gene tree in one common coordinate system so that
trees over different species subsets are comparable.

`readTrees` picks a master tree containing all species and enumerates every
**path** in it — every ordered pair `(descendant, ancestor)` where the second
node lies on the route from the first to the root:

```r
ap <- allPathsTT(master)     # ap$dist, ap$nodeId, ap$matIndex
```

`ap$matIndex[d, a]` gives the column of the `paths` matrix that holds path
`d -> a`, and is `NA` where `a` is not an ancestor of `d`. For the 62-taxon
mammal master this is 800 columns.

Each gene tree then contributes one row. A gene tree missing species is a
*contraction* of the master: some master nodes survive, others collapse. Paths
between two surviving nodes are computable; paths with a collapsed endpoint are
genuinely missing and are left `NA`. That is the intended and only legitimate
source of `NA`.

Mapping a gene tree's nodes onto master node numbers is done by
`matchAllnodesTT` (`R/RERfuncs.R:3749`):

```r
matchAllnodesTT <- function (tree, masterTree) {
  index = KeptVerts(masterTree, TipLabels(masterTree) %in% tree$tip.label)
  key   = which(index)
  cbind(seq_along(key), key)
}
```

**This mapping is purely positional.** It asserts that gene-tree vertex `i` is
the `i`-th surviving master vertex. Nothing about the gene tree is inspected
except its tip labels.

## 2. The bug

The positional assertion holds only if the gene tree's vertex numbering equals
the master's numbering restricted to the surviving vertices. That requires
three things:

1. the same tip order — `RenumberTips`
2. preorder traversal — `Preorder`
3. **the same root** — because preorder numbering is defined relative to a root

The released `readTrees` (`daadf41`) did the first two and not the third:

```r
trees[[i/2]] = unroot(read.tree(text = tmp[i]))    # per gene
...
master = Preorder(master)                          # master left unrooted
for (i in 1:treesObj$numTrees) {
  treesObj$trees[[i]] = RenumberTips(treesObj$trees[[i]], master$tip.label)
  treesObj$trees[[i]] = Preorder(treesObj$trees[[i]])
}
```

Everything is unrooted, so there is no canonical root to align to. Where `ape`
puts the basal trifurcation depends on how the newick string happened to be
written. Two files describing the *same* unrooted tree can number its internal
nodes differently, and the positional map then pairs up non-homologous nodes.

### Minimal demonstration

One unrooted 6-taxon tree, written two ways. Same topology (RF = 0), same
branch lengths (identical sorted edge-length vectors):

```
newick 1:  (((A:1,B:1):1,(C:1,D:1):1):1,(E:1,F:1):1);
newick 2:  (((C:1,D:1):1,(E:1,F:1):2):1,A:1,B:1);
```

All six species are present in both, so `KeptVerts` retains every vertex and
the map is the identity — gene node `k` is asserted to be master node `k`.
Under the old preparation (`unroot` + `RenumberTips` + `Preorder`):

```
node  7 :  master = ABCDEF   gene = ABCDEF   ok
node  8 :  master = AB       gene = CDEF     <-- MISMATCH
node  9 :  master = CD       gene = CD       ok
node 10 :  master = EF       gene = EF       ok
```

Master node 8 is the clade `AB`; the gene tree's node 8 is `CDEF`. The path
lengths written into the `AB` columns are measured from an unrelated part of
the tree.

Under the new preparation (`prepareTreeForTT` = `rootLikeMaster` +
`RenumberTips` + `Preorder`) against a rooted master, all five internal nodes
match.

## 3. What the failure looks like in practice

Three distinct outcomes, in increasing order of severity.

### (a) Computable paths dropped to NA

When the mismapped node pair happens not to be an ancestor pair in the master,
`matIndex` returns `NA` and the value is discarded even though the path is
perfectly computable. On a synthetic 12-taxon master with 300 random species
subsets, gene trees arriving with an independent root placement:

| | old | new |
|---|---|---|
| computable paths per gene | 51.9 | 66.5 |
| paths actually filled | 4.9 | 66.5 |
| computable paths lost to NA | **47.1 (91%)** | **0** |
| gene trees losing ≥1 path | 271 / 300 | 0 / 300 |

### (b) False "discordant tree topology" rejection

`allPathsMasterRelativeTT` treated any `NA` in the lookup as evidence that the
gene tree's topology disagrees with the master. The topology agrees exactly —
only the newick's root placement differs — so the diagnosis is wrong and the
whole gene is thrown away.

### (c) Rejected genes silently became rows of zeros

The rejection path is worse than the message suggests (`daadf41`,
`R/RERfuncs.R`):

```r
vals = double(length(masterTreePaths$dist))     # all ZEROS
if (sum(is.na(ii)) > 0 & !is.null(i)) {
  message("warning: discordant tree topology in tree ", i, ", returning NA row")
  return(vals)                                  # returns zeros, not NA
}
vals[] = NA                                     # only reached if NOT rejected
vals[ii] = treePaths$dist
```

`double(n)` is a vector of zeros. `vals[] = NA` sits *after* the early return,
so a rejected gene gets a row of **zeros, not `NA`s**. Those zeros pass every
`!is.na()` filter downstream and enter RER estimation as genuine zero-length
branches.

## 4. Scale on real data

The bundled `ext/subsetMammalGeneTrees.txt` does **not** trigger the bug, and
this is the normal case rather than a lucky one: all 500 trees came out of a
single pipeline and were written with a consistent root placement, so the
positional map happens to be correct and the old code loses nothing:

| | old | new |
|---|---|---|
| computable paths per gene | 1188.0 | 689.1 |
| lost to NA | 0 | 0 |
| spurious values | 0 | 0 |

The bug is latent, not absent. Re-rooting each of those same 500 trees at a
random tip — which changes only where the basal fork sits, leaving topology and
every branch length untouched (max RF = 0, identical sorted edge lengths) — and
feeding the rewritten file back in:

| | old | new |
|---|---|---|
| gene rows that became **all zeros** | **485 / 500** | 0 |
| spurious values per gene | 216.4 | 0 |
| share of filled values that are bogus | 15.8% | 0% |
| paths matrix dimensions | 500×1401 → 500×**1377** | 500×800 → 500×**800** |
| NA pattern vs. original run | changed | **identical** |

The old pipeline's *column space itself* changed (1401 → 1377 paths) because
the master tree is just the first complete gene tree, so re-rooting it changed
the set of enumerated paths. The new pipeline reproduces the master edge matrix,
the tip order, the column count and the NA pattern exactly.

## 5. The fix

```r
prepareTreeForTT <- function(tree, master) {
  tree <- rootLikeMaster(tree, master)          # NEW: align root position
  tree <- RenumberTips(tree, master$tip.label)
  Preorder(tree)
}
```

`rootLikeMaster` roots the gene tree on the edge separating the taxa that fall
on one side of the master's root, restricted to the species the gene tree
actually has. With the master kept rooted and every gene tree rooted
homologously, preorder numbering is now well defined and the positional map in
`matchAllnodesTT` is valid by construction.

`hasConcordantTopology` replaces the `NA`-count heuristic with an actual
topology comparison (prune both to shared taxa, canonicalise, compare newick),
so genuinely discordant trees are still caught — and now correctly return `NA`,
not zeros.

## 6. The root-edge split

The root is a point inserted *on* an edge, so how that edge's length is divided
between the root's two child branches is not determined by the unrooted tree.
`TreeTools::RootTree` and `phangorn::midpoint` both resolve this by putting the
whole length on one side and **zero** on the other:

```
input (unrooted):  ((A:1,B:1):1,(C:1,D:1):1,(E:1,F:1):2);
after RootTree  :  (((A:1,B:1):1,(C:1,D:1):1):2,(E:1,F:1):0);
```

Total length is conserved, so every pairwise distance is correct and nothing
looks malformed. But two things are wrong with it here:

1. A zero-length branch is a branch on which no change can occur — `expm(Q * 0)`
   is the identity — which misspecifies any model fitted on the rooted tree.
   `categoricalPermulations` failed outright on this: transition probability
   exactly 0 across that branch drove the annealing ratio to `0/0`.
2. **Which** side receives the zero depends on where the input newick happened
   to be rooted — reintroducing precisely the representation-dependence this
   package roots trees to eliminate (section 2).

It was not an edge case. Every gene tree came out of `rootLikeMaster` with a
zero-length root child, at the same edge index, so the master — whose branch
lengths are the mean of the complete gene trees' — inherited it deterministically.

`balanceRootEdges()` splits the edge evenly instead, placing the root at its
midpoint. The total is unchanged, so distances between all other nodes are
preserved, and the split is canonical rather than input-dependent. It is applied
to gene trees and the master (distance-like edge lengths), and deliberately
**not** to phenotype or trait trees, whose edge lengths carry trait values that
must not be averaged across the two halves of the split branch.

With this, path construction is exactly invariant to newick root placement:

| | `0 / full` split | even split |
|---|---|---|
| columns differing under re-rooting | 150 (122 root-incident) | **0** |

## 7. Reproducing

```sh
Rscript tests/reproduce_rooting_bug.R           # mechanism + new-pipeline invariance
```

For the old-vs-new comparison, build the pre-branch code in a worktree:

```sh
git worktree add /tmp/rer-old daadf41
R CMD INSTALL -l /tmp/rer-oldlib /tmp/rer-old
```

then run `readTrees` from each library on the original and re-rooted tree files.
