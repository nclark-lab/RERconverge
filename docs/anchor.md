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

* **Master.** Rooted on the edge from the anchor to its largest side `A`, with
  the root placed at the anchor node: a zero-length stem. The root stays binary.
  `TreeTools::KeepTip()` mis-sums edge lengths when pruning a side of a
  multifurcating root, so a multifurcating root is not used.
* **Gene with species on every side of the anchor.** It is rooted on the same
  edge, with the whole edge on side `A` (`anchorRootSplit()`). Its root coincides
  with the anchor node, so every filled column is a real distance.
* **Gene missing a side.** It is rooted on the master root restricted to its
  species, with the root edge split evenly (`balanceRootEdges()`). Only the
  columns ending at that root depend on the split, and their sum is the gene's
  edge length.
* **`anchor = "root"`** keeps the root of a rooted `masterTree`, as on
  `feature-rooted-trees`.

The anchor is stored in `treesObj$anchor` and on the master
(`masterTree$anchorSides`). Trait trees prepared with `prepareTreeForTT()` use
the same root placement.

## Discordance is decided by topology, not by lookup

`treeTopologyStatus()` compares each gene's bipartitions with the master pruned
to its species, using TreeTools' split code. Each gene gets a status:

* `"ok"` — same unrooted topology;
* `"unresolved"` — a polytomy in the gene;
* `"discordant"` — a split the master does not have;
* `"species_not_in_master"`, `"duplicate_tips"`.

Only `"ok"` genes are mapped; the rest get an all-`NA` row and are counted in a
message. Statuses are in `treesObj$treeStatus`.

For an `"ok"` gene, every path must have a master column. A missing column means
the tree was not rooted and numbered like the master, and
`allPathsMasterRelativeTT()` stops with an internal error. It never reports
discordance.

## Small trees

Gene trees with fewer than `minTreeSpecies` species (default 10, after
`useSpecies`) are dropped at input, with a message. They have no power for RER
analysis, and dropping them avoids special handling of 3- and 4-species trees.

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

## Tests

`tests/testthat/test_anchor.R`, with ape-only truth in
`tests/testthat/helper-anchor.R`:

* **Schema.** The master is rooted at the anchor with a zero stem; the columns
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
