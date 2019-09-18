#Modify initial path and inst/extdata on install (from where do tests run?)
#To do: Find a way to suppress printed output but not
repodir = "~/repos/RERconverge"
source(paste(repodir,"R/estimateTreeFuncs.R",sep="/"), echo=TRUE)
taln = paste(repodir,"inst/extdata/eg_alns/A1BG.phy", sep="/")
ttree = paste(repodir,"inst/extdata/MarineTreeBin.txt", sep="/")
talndir = paste(repodir,"inst/extdata/eg_alns", sep="/")
toutfile = paste(repodir,"inst/extdata/tmp_esttrees.txt", sep="/")
#test one tree first
tptree = estimatePhangornTree(taln, ttree, submodel="JC", type="DNA", k=1)
#then test all of the trees (reduce this to 10 trees)
tptrees = estimatePhangornTreeAll(alndir=talndir, pattern="*.phy", treefile=ttree,
                                  output.file=toutfile,
                                  submodel="JC", type = "DNA", format = "fasta", k=1)
