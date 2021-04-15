
library(stringr)

#read filenames
fof=read.table("../HairlessBSmodels/alnfof.txt", stringsAsFactors = F)
fof=fof$V1
fof=unlist(lapply(fof, function(x){paste0(treedir,strsplit(x, split="[.]")[[1]][1], ".tre")}))

#location of trimmed trees
treedir="../HairlessBSmodels/trimmedtrees/"

#desired location of foreground-specified trees
numdir="../HairlessBSmodels/Numberedtrees/"

#master tree
master=read.tree(text=readLines("../HairlessBSmodels/meredithplustree62sphg19rooted"))

#foreground species
fg=c("dasNov3", "triMan1", "loxAfr3", "cerSim1", "susScr3", "turTru2", "orcOrc1", "odoRosDi", "hetGla2", "hg19")
#check for typos - these should give the same number
# length(fg)
# sum(fg %in% master$tip.label)



####################################

#add numbers
for(f in fof){
  tree=readLines(f)
  for(s in fg){
    tree=str_replace(tree, s, paste0(s, " #1", collapse=""))
  }
  outfile=paste0(numdir,strsplit(f, "[/]")[[1]][4], collapse="")
  write.table(tree, file=outfile, quote = F, col.names = F, row.names = F)
}






