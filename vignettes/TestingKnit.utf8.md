---
title: "TestingKnit"
author: "Wynn Meyer"
date: "4/3/2019"
output: pdf_document
---



## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:


```r
summary(cars)
```

```
##      speed           dist       
##  Min.   : 4.0   Min.   :  2.00  
##  1st Qu.:12.0   1st Qu.: 26.00  
##  Median :15.0   Median : 36.00  
##  Mean   :15.4   Mean   : 42.98  
##  3rd Qu.:19.0   3rd Qu.: 56.00  
##  Max.   :25.0   Max.   :120.00
```

## Including Plots

You can also embed plots, for example:

![](/private/var/folders/t3/c42vfq915ljf65t3gslyfzzm0000gp/T/RtmppT9dtm/preview-3ff721f5f5c.dir/TestingKnit_files/figure-latex/pressure-1.pdf)<!-- --> 

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.


## My own work

Attempting to load and plot RER data.


```r
library("RERconverge")
data("toyTrees") 
data("mamRERw")
#repopath = "~/repos"
#source(paste(repopath,"/R/RERfuncs.R",sep=""))
#source(paste(repopath,"/R/plottingFuncs.R",sep="")) #for troubleshooting margins
```

Here is where I have issues with figure margins.


```r
#plot RERs as tree
#par(mfrow=c(1,1))
newbend3rers = treePlotRers(treesObj=toyTrees, rermat=mamRERw, index="BEND3", 
                            type="c", nlevels=9, figwid=10)
```

![](/private/var/folders/t3/c42vfq915ljf65t3gslyfzzm0000gp/T/RtmppT9dtm/preview-3ff721f5f5c.dir/TestingKnit_files/figure-latex/fig2c-1.pdf)<!-- --> 

```r
#treePlotRers(treesObj=toyTrees, rermat=mamRERw, index="BEND3", type="c", nlevels=9)
```
