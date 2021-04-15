
#get filenames
fns=list.files("../HairlessBSmodels/Checkfortreewidepositiveselection/outM1round2/")
m1fns=paste0("../HairlessBSmodels/Checkfortreewidepositiveselection/outM1round2/", fns)
m2fns=paste0("../HairlessBSmodels/Checkfortreewidepositiveselection/outM2round2/", fns)
m8afns=paste0("../HairlessBSmodels/Checkfortreewidepositiveselection/outM8Around2/", fns)
m8fns=paste0("../HairlessBSmodels/Checkfortreewidepositiveselection/outM8round2/", fns)

#output file
outfile="../HairlessBSmodels/BSresultsRound2.csv"

#extract lnL statistics

allm1=c()
allm2=c()
allm8a=c()
allm8=c()

count=1
while(count<=length(m1fns)){
  curm1=readLines(m1fns[count])
  curm2=readLines(m2fns[count])
  curm8a=readLines(m8afns[count])
  curm8=readLines(m8fns[count])
  
  m1=strsplit(curm1[which(grepl("lnl", curm1, ignore.case = T))], split="\\s+")[[1]][5]
  m2=strsplit(curm2[which(grepl("lnl", curm2, ignore.case = T))], split="\\s+")[[1]][5]
  m8a=strsplit(curm8a[which(grepl("lnl", curm8a, ignore.case = T))], split="\\s+")[[1]][5]
  m8=strsplit(curm8[which(grepl("lnl", curm8, ignore.case = T))], split="\\s+")[[1]][5]
  
  allm1=c(allm1, m1)
  allm2=c(allm2, m2)
  allm8a=c(allm8a, m8a)
  allm8=c(allm8, m8)
  count=count+1
}


#calculate log-likelihood ratio
#M2-M1 #simpler discrete model
#M8-M8A #more complex continuous model
#p-value: chi-square, right-tail, 1 df

bsdf=data.frame(m1=as.numeric(allm1), m2=as.numeric(allm2), m8=as.numeric(allm8), m8a=as.numeric(allm8a))
rownames(bsdf)=fns

bsdf$LRTm2m1=bsdf$m2-bsdf$m1
bsdf$LRTm8m8a=bsdf$m8-bsdf$m8a

bsdf$pvalm2m1=pchisq(2*bsdf$LRTm2m1, df=1, lower.tail = F)
bsdf$pvalm8m8a=pchisq(2*bsdf$LRTm8m8a, df=1, lower.tail = F)
write.csv(bsdf, outfile)



