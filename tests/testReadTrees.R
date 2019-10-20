#When running tests, revise the paths to the data.
#Run testEstTrees first (ensures tmp_esttrees is created)
source("tests/testEstTrees.R")
tmptrees = "inst/extdata/tmp_esttrees.txt"
masttree = "inst/extdata/MarineTreeBin.txt"
mtree = read.tree(masttree)
mtree$edge.length = runif(length(mtree$edge.length))
source("R/RERfuncs.R")
trtrees = readTrees(tmptrees, masterTree=mtree, reorient=T,outgroup="ornAna1")
