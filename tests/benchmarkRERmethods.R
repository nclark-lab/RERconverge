# Benchmark methods to calculate RERs against original and updated methods of RERconverge on 
# correlations with subterranean phenotype for eye-specific genes


if(!require(PRROC)){
     install.packages('PRROC')
     library(PRROC)
}

source('R/RERfuncs.R')

rerpath = find.package('RERconverge')
rerpath <- '.'
originalupdatedprsfile = "OriginalUpdatedRERsMoleAccelerationPrecisionRecall.rds"
originalupdatedprs=readRDS(paste(rerpath,"/ext/",originalupdatedprsfile,sep=""))
eyemastergenesfile = "eyespecificgenes.master.list"
eyemastergenes=readLines(paste(rerpath,"/ext/",eyemastergenesfile,sep=""))
moles <- c('Naked_mole_rat','Star_nosed_mole','Blind_mole_rat','Cape_golden_mole')

#user input begin
pathtoreesobjectfile <- '~/Documents/rermethods/data/mamm63nt.trees.rds'
pathtoRERmatrixobject <- '~/Documents/rermethods/data/mamm63nt.trees.plac.scaledrers.sqrt.wt.rds'
foregroundforbenchmark <- moles[-4] #change to moles if all moles are to be used
treesObj <- readRDS(pathtoreesobjectfile)
rermat <- readRDS(pathtoRERmatrixobject)
figureoutputfilename <- 'benchmark_moleAcceleration_eyemastergenes.pdf'
#user input end

#Run rest of the code to generate the Precision-recall plot

calcprc <- function(scores,labels){
     fg <- scores[labels == 1]
     bg <- scores[labels == 0]
     pr <- pr.curve(scores.class0 = fg[!is.na(fg)], scores.class1 = bg[!is.na(bg)], curve = T)
     return(pr)
}
addontoprecrecplot <- function(rermat, foreground, treesObj, positive.gene.list){
     phenvfgd <- foreground2Paths(foreground, treesObj,  clade = 'terminal')
     rermat.fgdcorr <- correlateWithBinaryPhenotype(rermat, phenvfgd) %>% 
          mutate(gene = rownames(rermat)) %>%
          rowwise() %>%
          mutate(score = -1*sign(Rho)*log(P)) %>% arrange(desc(score))
     rermat.fgdcorr.prs <- calcprc(rermat.fgdcorr$score, 1*(rermat.fgdcorr$gene %in% positive.gene.list))
     lines(rermat.fgdcorr.prs$curve[,1],
           rermat.fgdcorr.prs$curve[,2], lty = 1, lwd = 2, col = 'blue')
     legend('top',c('Your method'),col=c('Blue'),pch=19,title='',bty='n')
}


org.mole.prs <- originalupdatedprs$original.mole
org.control.prs <- originalupdatedprs$original.control
upd.mole.prs <- originalupdatedprs$updated.mole
upd.control.prs <- originalupdatedprs$updated.control
pdf(paste0('tests/',figureoutputfilename), width = 5, height = 5)
plot(org.mole.prs$curve[,1], org.mole.prs$curve[,2], type = 'l', lty = 1,xlim = c(0,0.5),lwd=2,
     xlab = 'Recall',ylab = 'Precision')
title('Eye-specific master genes', adj = 0)
lines(org.control.prs$curve[,1], org.control.prs$curve[,2], lty = 2,lwd=2)
lines(upd.control.prs$curve[,1], upd.control.prs$curve[,2], lty = 2, col = 'red',lwd=2)
lines(upd.mole.prs$curve[,1], upd.mole.prs$curve[,2], lty = 1, col = 'red',lwd=2)
legend('topright',c('Mole','Control'),lty=c(1,2),title='Foreground', bty = 'n')
legend('right',c('Original','Updated'),col=c('Black','Red'),pch=19,title='RERconverge',bty='n')
addontoprecrecplot(rermat = rermat, 
                   foreground = foregroundforbenchmark, 
                   treesObj = treesObj, 
                   positive.gene.list = eyemastergenes)
dev.off()