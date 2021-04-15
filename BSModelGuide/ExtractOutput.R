
#get filenames
fns=list.files("../HairlessBSmodels/outBS1")
bs1fns=paste0("../HairlessBSmodels/outBS1/", fns)
bs2fns=paste0("../HairlessBSmodels/outBS2/", fns)
m1fns=paste0("../HairlessBSmodels/outM1/", fns)

#output file
outfile="../HairlessBSmodels/BSresults.csv"

#extract lnL statistics

allbs1=c()
allbs2=c()
allm1=c()
count=1
while(count<=length(bs1fns)){
  curbs1=readLines(bs1fns[count])
  curbs2=readLines(bs2fns[count])
  curm1=readLines(m1fns[count])
  bs1=strsplit(curbs1[which(grepl("lnl", curbs1, ignore.case = T))], split="\\s+")[[1]][5]
  bs2=strsplit(curbs2[which(grepl("lnl", curbs2, ignore.case = T))], split="\\s+")[[1]][5]
  m1=strsplit(curm1[which(grepl("lnl", curm1, ignore.case = T))], split="\\s+")[[1]][5]
  allbs1=c(allbs1, bs1)
  allbs2=c(allbs2, bs2)
  allm1=c(allm1, m1)
  count=count+1
}


#calculate log-likelihood ratio
#BS1-M1 #foreground acceleration
#BS2-BS1 #foreground positive selection
#p-value: chi-square, right-tail, 1 df

bsdf=data.frame(m1=as.numeric(allm1), bs1=as.numeric(allbs1), bs2=as.numeric(allbs2))
rownames(bsdf)=fns

bsdf$LRTbs1m1=bsdf$bs1-bsdf$m1
bsdf$LRTbs2bs1=bsdf$bs2-bsdf$bs1

bsdf$AccPval=pchisq(2*bsdf$LRTbs1m1, df=1, lower.tail = F)
bsdf$PosSelPval=pchisq(2*bsdf$LRTbs2bs1, df=1, lower.tail = F)
write.csv(bsdf, outfile)



