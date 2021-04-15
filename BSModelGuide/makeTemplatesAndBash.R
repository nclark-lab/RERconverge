
bashfn="../HairlessBSmodels/runallBS.bash"

template1=readLines("../HairlessBSmodels/Round1templates/codeml_templateBS1")
template2=readLines("../HairlessBSmodels/Round1templates/codeml_templateBS2")
template3=readLines("../HairlessBSmodels/Round1templates/codeml_templateM1")

#read filenames
fof=read.table("../HairlessBSmodels/alnfof.txt", stringsAsFactors = F)
fof=fof$V1
#specify locations RELATIVE TO BASH FILE
foftree=unlist(lapply(fof, function(x){paste0("./Numberedtrees/",strsplit(x, split="[.]")[[1]][1], ".tre")}))
fofaln=unlist(lapply(fof, function(x){paste0("./trimmedalns/",x)}))
fofsimple=unlist(lapply(fof, function(x){strsplit(x, "[.]")[[1]][1]}))

#filled in templates location
templatedir="../HairlessBSmodels/writtentemplates/"
#template location relative to bash file
reltemplatedir="./writtentemplates/"

count=1
while(count<=length(foftree)){
  tf=foftree[count]
  af=fofaln[count]
  sf=fofsimple[count]
  
  template1[1]=paste0("      seqfile =  ",af,"  * sequence data filename", collapse="")
  template2[1]=paste0("      seqfile =  ",af,"  * sequence data filename", collapse="")
  template3[1]=paste0("      seqfile =  ",af,"  * sequence data filename", collapse="")
  
  template1[2]=paste0("     treefile =  ",tf,"     * tree structure file name", collapse="")
  template2[2]=paste0("     treefile =  ",tf,"     * tree structure file name", collapse="")
  template3[2]=paste0("     treefile =  ",tf,"     * tree structure file name", collapse="")
  
  template1[3]=paste0("      outfile =  ",paste0("./outBS1/",sf, collapse=""),"       * main result file name", collapse="")
  template2[3]=paste0("      outfile =  ",paste0("./outBS2/",sf, collapse=""),"       * main result file name", collapse="")
  template3[3]=paste0("      outfile =  ",paste0("./outM1/",sf, collapse=""),"       * main result file name", collapse="")
  
  t1=paste0(templatedir, "BS1", sf)
  t2=paste0(templatedir, "BS2", sf)
  t3=paste0(templatedir, "M1", sf)
  
  conn=file(t1)
  writeLines(template1, conn)
  close(conn)
  
  conn=file(t2)
  writeLines(template2, conn)
  close(conn)
  
  conn=file(t3)
  writeLines(template3, conn)
  close(conn)
  
  count=count+1
}



#make bash file to run all
allfiles=list.files(templatedir)
allfiles=unlist(lapply(allfiles, function(x){paste0(reltemplatedir, x, collapse="")}))

bashlines=c()
for(f in allfiles){
  st=paste0("nohup codeml ",f," &", collapse="")
  bashlines=c(bashlines,st)
}

conn=file(bashfn)
writeLines(bashlines, conn)
close(conn)








