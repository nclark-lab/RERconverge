
bashfn="../HairlessBSmodels/Checkfortreewidepositiveselection/round2bash.sh"

#alns and trees
alnfns=list.files("../HairlessBSmodels/Checkfortreewidepositiveselection/Round2trimmedalns/")
treefns=list.files("../HairlessBSmodels/Checkfortreewidepositiveselection/Round2trimmedtrees/")

template1=readLines("../HairlessBSmodels/Checkfortreewidepositiveselection/Round2templates/codeml_templateM1")
template2=readLines("../HairlessBSmodels/Checkfortreewidepositiveselection/Round2templates/codeml_templateM2")
template3=readLines("../HairlessBSmodels/Checkfortreewidepositiveselection/Round2templates/codeml_templateM8")
template4=readLines("../HairlessBSmodels/Checkfortreewidepositiveselection/Round2templates/codeml_templateM8A")

#location to write template files
templatedir="../HairlessBSmodels/Checkfortreewidepositiveselection/writtenround2templates/"

#written template location relative to bash file
reltemplatedir="./writtenround2templates/"
#tree and aln locations relative to bash file
af="./Round2trimmedalns/"
tf="./Round2trimmedtrees/"

count=1
while(count<=length(alnfns)){
  n=strsplit(alnfns[count], split="[.]")[[1]][1]
  
  template1[1]=paste0("      seqfile =  ", paste0(af, n, ".phy", collapse=""), "  * sequence data filename")
  template1[2]=paste0("     treefile =  ",paste0(tf, n, ".tre") ,"     * tree structure file name")
  template1[3]=paste0("      outfile =  ",paste0("./outM1round2/", n),"       * main result file name")
  
  template2[1]=paste0("      seqfile =  ", paste0(af, n, ".phy", collapse=""), "  * sequence data filename")
  template2[2]=paste0("     treefile =  ",paste0(tf, n, ".tre") ,"     * tree structure file name")
  template2[3]=paste0("      outfile =  ",paste0("./outM2round2/", n),"       * main result file name")
  
  template3[1]=paste0("      seqfile =  ", paste0(af, n, ".phy", collapse=""), "  * sequence data filename")
  template3[2]=paste0("     treefile =  ",paste0(tf, n, ".tre") ,"     * tree structure file name")
  template3[3]=paste0("      outfile =  ",paste0("./outM8round2/", n),"       * main result file name")
  
  template4[1]=paste0("      seqfile =  ", paste0(af, n, ".phy", collapse=""), "  * sequence data filename")
  template4[2]=paste0("     treefile =  ",paste0(tf, n, ".tre") ,"     * tree structure file name")
  template4[3]=paste0("      outfile =  ",paste0("./outM8Around2/", n),"       * main result file name")
  
  t1=paste0(templatedir, "M1", n)
  t2=paste0(templatedir, "M2", n)
  t3=paste0(templatedir, "M8", n)
  t4=paste0(templatedir, "M8A", n)
  
  conn=file(t1)
  writeLines(template1, conn)
  close(conn)
  
  conn=file(t2)
  writeLines(template2, conn)
  close(conn)
  
  conn=file(t3)
  writeLines(template3, conn)
  close(conn)
  
  conn=file(t4)
  writeLines(template4, conn)
  close(conn)
  
  count=count+1
}




#make bash file to run all
allfiles=list.files(templatedir)
allfiles=unlist(lapply(allfiles, function(x){paste0(reltemplatedir, x, collapse="")}))

bashlines=c()
count=1
for(f in allfiles){
  if(count %% 4==0){
    st=paste0("nohup codeml ",f, collapse="")
  }else{
    st=paste0("nohup codeml ",f," &", collapse="")
  }
  bashlines=c(bashlines,st)
  count=count+1
}

conn=file(bashfn)
writeLines(bashlines, conn)
close(conn)






