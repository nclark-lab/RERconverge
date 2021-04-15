library(phytools)
library(phangorn)
library(seqinr)

#read filenames
fof=read.table("../HairlessBSmodels/alnfof.txt", stringsAsFactors = F)
fof=fof$V1

#location of alignments
alndir="../HairlessBSmodels/mamm63nt.alns/trimallcds_wSpalax/"

#master tree
master=read.tree(text=readLines("../HairlessBSmodels/meredithplustree62sphg19rooted"))

#desired directory for trimmed alignments
finaldir="../HairlessBSmodels/trimmedalns/"

#desired directory for trimmed trees
treedir="../HairlessBSmodels/trimmedtrees/"

#all species to include
specstouse=c("dasNov3", "triMan1", "loxAfr3", "oryAfe1", "eleEdw1",
             "sorAra2", "conCri1", "pteAle1", "myoDav1", "equCab2",
             "cerSim1", "vicPac2", "susScr3", "turTru2", "orcOrc1",
             "oviAri3", "bosTau7", "felCat5", "musFur1", "odoRosDi",
             "lepWed1", "tupChi1", "oryCun2", "jacJac1", "mm10",
             "criGri1", "speTri2", "hetGla2", "cavPor3", "otoGar3",
             "calJac3", "rheMac3", "chlSab1", "panTro4", "hg19")
#check for typos - these should give the same number
# length(specstouse)
# sum(specstouse %in% master$tip.label)



################################################################################



#make trimmed aligns and pruned trees
for(f in fof){
  #align
  curf=paste0(alndir, f, collapse="")
  curaln=read.fasta(curf, as.string = T)
  l=nchar(curaln[[1]][1])+100
  curaln=curaln[names(curaln) %in% specstouse]
  outfile=paste0(finaldir, f)
  write.fasta(sequences = curaln, names = names(curaln),file.out = outfile, as.string = T, nbchar = l)
  #tree
  specsfortree=names(curaln)
  totrim=master$tip.label[which(!(master$tip.label %in% specsfortree))]
  newtree=master
  newtree=drop.tip(master, totrim)
  n=paste0(strsplit(f, split="[.]")[[1]][1], ".tre", collapse="")
  treeoutfile=paste0(treedir, n)
  write.tree(newtree, file=treeoutfile)
}










