# Anchored path representation

This branch keeps the `paths` matrix: one column per ancestor–descendant pair of
the master tree. It changes where the master is rooted, how gene trees are
rooted against it, and how discordance is detected.

## Why an anchor

A gene tree is unrooted. For a column to hold a real distance, both of its
endpoints must be nodes of the gene's unrooted tree. The master's root is such a
node only when it is a real internal node **and** the gene has species on every
side of it. Otherwise the gene's root is a point on an edge, and every path
ending there depends on an arbitrary split of that edge.

* Rooting the master on an edge (the previous `feature-rooted-trees` design)
  makes every gene's root a point on an edge. Every gene has arbitrary root
  values.
* The original unrooted design anchored the columns at the basal node that
  `unroot()` happened to leave. It had no node to store the crossing edge when a
  gene missed a side, and it mapped nodes by position, which depends on how the
  newick was written.

`readTrees(anchor = "auto")` picks the internal node with the most genes covering
all of its sides. On the 466-species, 16,209-gene issue dataset, the best node
leaves 15 genes (0.1%) without every side. The monotreme trifurcation that
`unroot()` leaves would leave 8,440 (52%).

## Representation

* **Master.** Rooted at the anchor node itself, as a multifurcating root with
  one child per side (`anchorMaster()`, `rootAtNodeJoining()`). No edge or
  vertex is added, so there is no zero-length stem, no all-zero column and no
  duplicated columns.
* **Gene with species on every side of the anchor.** It is rooted at the same
  node (`rootLikeMaster()`). Every filled column is a real distance.
* **Gene missing a side.** The master's root, pruned to the gene's species, is
  then a point on an edge. The gene is rooted on that edge and the edge is split
  evenly (`balanceRootEdges()`). Only the columns ending at that point depend on
  the split, and their sum is the gene's edge length.
* **`anchor = "root"`** keeps the root of a rooted `masterTree`, as on
  `feature-rooted-trees`.
* **Tie-break.** The anchor is the node with the most genes covering all of its
  sides. Ties are broken by side sizes and a key built from the unordered sides,
  so the choice does not depend on how the master newick is written.

TreeTools' pruning helpers are not used on the master.
`TreeTools::KeepTip()` and `TreeTools::KeptVerts()` suppress the wrong node,
and `KeepTip()` also mis-sums edge lengths, when a side of a multifurcating
root is dropped. Instead:

* `rootLikeMaster()` and `treeTopologyStatus()` prune with `ape::keep.tip()`;
  it is also about ten times faster on a 466-species master.
* `matchAllnodesTT()` computes the kept master vertices itself. A tip is kept if
  present; an internal node is kept if at least two of its child subtrees
  contain the tree's species.

The anchor is stored in `treesObj$anchor`. Trait trees prepared with
`prepareTreeForTT()` are rooted the same way and their edge values are not
changed.

## Discordance is decided by topology, not by lookup

`treeTopologyStatus()` compares each gene's bipartitions with the master pruned
to its species. Each gene gets a status:

* `"ok"` — same unrooted topology;
* `"unresolved"` — a polytomy in the gene;
* `"discordant"` — a split the master does not have;
* `"species_not_in_master"`, `"duplicate_tips"`.

Genes that are not `"ok"` are dropped at input, like trees that are too small,
and counted in a message. `treesObj$dropped` records each dropped gene and its
reason (`"few_species"` or its status), and `treesObj$treeStatus` is `"ok"`
throughout. Nothing that was dropped contributes to the species report, the
anchor choice, the master branch length estimate, or imputation and PC
normalisation.

Classification happens before the master is anchored, since the topology check
does not depend on rooting, so the anchor is chosen from the genes that are kept.

The RER export functions (`returnRersAsTree()` and friends) still tolerate a
`treesObj` that carries trees without paths, e.g. one saved before these trees
were dropped: they return those trees with `NA` edge lengths and warn.

For an `"ok"` gene, every path must have a master column. A missing column means
the tree was not rooted and numbered like the master, and
`allPathsMasterRelativeTT()` stops with an internal error. It never reports
discordance.

## Master branch lengths

Trait analysis reconstructs ancestral states on the master's branch lengths
(`char2Paths()` -> `edgeVars()` -> `phytools::fastAnc`), so by default they come
from the data: `masterBranchLengths = "estimate"` estimates each master edge
from the gene trees that have it, each gene scaled to its own total. `minSpecs`
optionally restricts which genes contribute, and fewer than `minTreesAll` usable
genes is an error rather than a master left without usable lengths.

Keeping an external tree's branch lengths is an explicit override,
`masterBranchLengths = "supplied"`, and that tree must have branch lengths. The
old `reestimateBranches` argument is deprecated: `FALSE` maps to `"supplied"`,
`TRUE` to `"estimate"`.

Master edges that are never a gene edge keep their supplied length, or 0 when
there is none.

## Small trees

Gene trees with fewer than `minTreeSpecies` species (default 10, after
`useSpecies`) are dropped at input, with a message. They have no power for RER
analysis, and dropping them avoids special handling of 3- and 4-species trees.

## Trait analysis

The anchor is a computational root: it is chosen from gene coverage and can move
when the gene set changes. Ancestral states, the direction of change and
"ancestral" clades are biological notions, so they must not depend on it.

`readTrees` therefore keeps a second rooting of the same master,
`treesObj$masterTreeRooted`, with the same (estimated) branch lengths:

* the rooting of the supplied `masterTree`, when there is one;
* otherwise the midpoint, with a warning.

`char2Paths()` reconstructs ancestral states on that rooting and writes the
per-branch values into the anchored columns:

* a branch's value is oriented away from the biological root, whichever way its
  anchored column runs;
* the one branch that contains that root has no single descendant end, so its
  value is the difference across the whole branch, oriented by the tip labels
  (away from the side holding the first label). Orienting it by the anchored
  column instead would make its sign depend on the anchor: an anchor on the other
  side of that branch traverses it the other way. It is the only value an anchor
  change could still alter, and on the 742-species set it did, flipping
  -0.0578 to +0.0578;
* nodes are matched between the two rootings by the bipartitions of their
  incident edges (`nodeIdentity()`, `mapNodesBetweenRootings()`), which does not
  depend on either rooting.

A master branch of length 0 -- identical sequences, or an edge no gene resolves
-- makes the reconstruction singular: a node whose children all sit at distance 0
weights their contrast by 1/(0 + 0), and the resulting NaN spreads to every state
in the tree, so the entire trait vector comes back NA with nothing said. Those
branches are given a length far below the smallest real one, with a warning, and
a non-finite reconstruction is an error rather than a silent all-NA result. On
the 742-species master this is not hypothetical: a 497-gene subset estimated
three zero-length cherries and voided all 19,645 trait values.

`tree2Paths()` follows the same rules as `readTrees` for a phenotype tree:

* topology is compared by bipartitions (`treeTopologyStatus()`), so a concordant
  tree written from another root is not rejected;
* a path with no master column is an internal error, not a silently dropped
  value;
* categorical values are the state of the branch below a node, so they are read
  before pruning and carried by branch, named by the tips below it. Where pruning
  merges branches the merged branch takes the most recent state, and the two
  branches at a tree's own root keep their own states. Continuous and binary
  values keep adding along merged branches.

## Other fixes on this branch

* `coreGetResiduals` prepares a pruned gene like the master before the positional
  edge lookup (it previously unrooted it, which broke internal-node indices). It
  skips genes without paths.
* `readTrees(minSpecs = )` reads each master edge from its own column. It
  previously used `ap$destinNode`, which does not exist.
* A master with species absent from every gene no longer breaks column naming
  (`namePathsWSpeciesTT`) or residuals.
* `rootLikeMaster` falls back to the complementary outgroup when
  `TreeTools::RootTree()` leaves a four-species tree unrooted.
* `concordant_trees` converts TreeTools-preordered trees before calling ape,
  which previously segfaulted.
* The legacy path consumers (`getProjectionPaths`, `correlateTreesAll`,
  `correlateTreesBinary`, `plotTreesBinary`, `plotContinuousChar`,
  `plotContinuousCharXY`) prepare their pruned trees instead of unrooting them,
  and take columns from `matIndex` by node pair instead of matching
  `namePaths()` against column names, which now hold species names.
  `getProjectionPaths` returned an all-`NA` matrix before this.
* The three `unroot(drop.tip(...))` calls are guarded with `apeOrder()`:
  `drop.tip` keeps TreeTools' preorder attribute when nothing is dropped, and
  `unroot` then produced an edge matrix with a self-loop and a tip with no edge.
* `correlateTreesBinary` compared its two column vectors with `all(ii=ii2)`, an
  assignment that is always `TRUE`.

## Cost

`getAllResiduals` computed a full weights matrix on every call -- a complete
unweighted regression over every gene, a lowess fit and diagnostic plots -- and
then discarded it whenever `use.weights = FALSE` or external weights were
supplied. It also built a full-size matrix of ones and squared it back through
`sqrt()`, materialised the inverse transform of the whole paths matrix only to
take column means of it, extracted every non-NA path value to take one quantile,
and copied the residual matrix twice more to normalise it.

None of that changes a number, so it is no longer paid for:

* weights are computed only when they are used
  (`transformPaths(computeWeights = )`);
* unit weights are carried as a flag and the slice each regression needs is built
  in the loop. Multiplying by `sqrt(1)` is exact, so the studentizing step is
  skipped rather than applied;
* `colStatBlocks()` takes column statistics a block of columns at a time, passing
  each column to the same function as before, so the values match
  `apply(mat, 2, f)` while no full-size temporary is built;
* `lowQuantileExact()` keeps only the smallest k values rather than all of them.
  It reproduces `quantile(type = 7)` exactly, including its
  `(1 - h) * lo + h * hi` form -- the algebraically equal `lo + h * (hi - lo)`
  rounds one ULP differently, and `cutoff` feeds a threshold comparison where a
  last bit can change which branches are masked;
* `getRMat()` applies the column statistics to the values it keeps instead of
  normalising a full copy of the matrix first.

On the 742-species set with 1500 genes this is 42.5 min and 6.3 GB before,
9.6 s and 3.4 GB after, and the results are bitwise identical: the RER matrices
and every column of the correlation table (`Rho`, `N`, `P`, `p.adj`). Checked
bitwise at 500 and 1500 genes unweighted, with computed weights, and with
external weights (which differ from the unweighted result, so the weighted path
is really exercised).

The weighted path itself is untouched and remains expensive -- it is dominated by
the lowess fit and the diagnostic plots inside `computeWeightsAllVar`, which cost
~13 min and ~5.9 GB on a 60 MB paths matrix. `impute = TRUE` and `n.pcs > 0` are
likewise unchanged, and were not part of this comparison.

## Tests

`tests/testthat/test_anchor.R`, with ape-only truth in
`tests/testthat/helper-anchor.R`:

* **Schema.** The master is rooted at the anchor node with no added edge; the columns
  are exactly the master's ancestor–descendant pairs; the anchor maximises gene
  coverage (brute force).
* **Master input.** The anchor and paths are identical whether the master is
  given rooted, unrooted or rewritten.
* **Values.** Every filled value equals the ape truth. Masters: random,
  caterpillar, balanced and small. Edge lengths: normal, 0/1e-8, 1e-12..10 and
  equal. Genes cover every side of the anchor or miss one.
* **Newick layout.** Rows are identical when the newick is rooted at a tip,
  rooted on an edge, unrooted or ladderized.
* **Pairs.** Branch vectors for two genes pruned to their common species, read
  from `paths`, equal direct pruning.
* **Pruning.** Values between real nodes do not move when species are dropped.
* **Status.** Discordant (tip swap, NNI) and unresolved genes are flagged with
  `NA` rows; concordant genes in any layout never are.
* **Lookup.** A failed lookup on a concordant tree is an error.
* **Inputs.** Small trees are dropped; master-only species; `minSpecs`
  re-estimation.
* **Residuals.** They use each gene's own branches, with and without
  `useSpecies`.
* **`anchor = "root"`.**
* **Regressions.** The four-species `RootTree` layout; `concordant_trees` on
  preordered trees.
