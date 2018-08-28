bdir = '../'
ctabfile = paste(bdir, 'SpeciesLists/SpeciesNameConversionTable0818.csv',sep='')
if (!file.exists(ctabfile)) {
  print(paste('Cant find species conversion table:',ctabfile,sep='/n'))
} else {
  cctab = read.csv(ctabfile, header=T, as.is=T) #does this make the df available in the workspace?
}
convertSpeciesNames <- function(tree1, tree2, ctab) {
  #convert species names in tree1 to the format of tree2, using ctab (a data frame)
  tree1conv = tree1
  for (t in 1:length(tree1$tip.label)) {
    if (tree1$tip.label[t] %in% tree2$tip.label == F) {
      #find the tip label somewhere in ctab
      ft = which(ctab==tree1$tip.label[t],arr.ind=T)
      if (dim(ft)[1] > 0){
        for (r in 1:length(unique(ft[,1]))) {
          rr = ctab[unique(ft[,1])[r],]
          if (sum(tree2$tip.label %in% rr) == 1) {
            tree1conv$tip.label[t] = tree2$tip.label[tree2$tip.label %in% rr]
          } else if (sum(tree2$tip.label %in% rr) > 1) {
            print(paste("Multiple options for",tree1$tip.label[t],"in tree2:"))
            print(tree2$tip.label[tree2$tip.label %in% rr])
          } else {
            print(paste(tree1$tip.label[t],"not found in tree2"))
          }
        }
      } else {
        print(paste(tree1$tip.label[t],"not found in conversion table"))
      }
    }
  }
  return(tree1conv)
}
