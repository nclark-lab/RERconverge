library(phytools)
library(phangorn)
library(seqinr)

#remove foreground species from trees
fg=c("dasNov3", "triMan1", "loxAfr3", "cerSim1", "susScr3", "turTru2", "orcOrc1", "odoRosDi", "hetGla2", "hg19")

#all species to include
specstouse=c("dasNov3", "triMan1", "loxAfr3", "oryAfe1", "eleEdw1",
             "sorAra2", "conCri1", "pteAle1", "myoDav1", "equCab2",
             "cerSim1", "vicPac2", "susScr3", "turTru2", "orcOrc1",
             "oviAri3", "bosTau7", "felCat5", "musFur1", "odoRosDi",
             "lepWed1", "tupChi1", "oryCun2", "jacJac1", "mm10",
             "criGri1", "speTri2", "hetGla2", "cavPor3", "otoGar3",
             "calJac3", "rheMac3", "chlSab1", "panTro4", "hg19")

#all non-foreground species
nonfg=setdiff(specstouse, fg)

#directory of trimmed alns round1
alnsdir="../HairlessBSmodels/trimmedalns/"
#directory of trimmed trees round1
treesdir="../HairlessBSmodels/trimmedtrees/"

#directory to put round2 trimmed alns
alnsoutdir="../HairlessBSmodels/Checkfortreewidepositiveselection/Round2trimmedalns/"
#directory to put round2 trimmed trees
treesoutdir="../HairlessBSmodels/Checkfortreewidepositiveselection/Round2trimmedtrees/"

#trim alns and trees
alnfns=list.files(alnsdir)
treefns=list.files(treesdir)

count=1
while(count<=alnfns){
  #aligns
  curf=paste0(alndir, alnfns[count], collapse="")
  curaln=read.fasta(curf, as.string = T)
  l=nchar(curaln[[1]][1])+100
  curaln=curaln[names(curaln) %in% nonfg]
  outfile=paste0(alnsoutdir, alnfns[count])
  write.fasta(sequences = curaln, names = names(curaln),file.out = outfile, as.string = T, nbchar = l)
  
  #trees
  newtree=read.tree(file=paste0(treedir, treefns[count], collapse=""))
  newtree=drop.tip(newtree, fg)
  treeoutfile=paste0(treesoutdir, treefns[count])
  write.tree(newtree, file=treeoutfile)
  
  count=count+1
}


#make templates for M1, M2, M8, M8A

#make bash file to run templates





