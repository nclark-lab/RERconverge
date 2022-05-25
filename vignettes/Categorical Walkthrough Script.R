# Walkthrough Script For Using New Categorical Functions

# Read in the inputs

# make sure you are in the directory that contains the files you are reading in

# read in the trees
trees=readRDS("zoonomiatrees.rds")

# read in the phenotype data
phenvals = readRDS("sleepPhenvals.rds")

# get the names of the species to include
allspecs = names(phenvals)

# read in the RERs from running getAllResiduals or run getAllResiduals
RERmat = readRDS("sleepPhenRERs.rds")
# or...
#RERmat = getAllResiduals(trees, useSpecies = allspecs)

# Generate paths using char2PathsCategorical (castor method, default)
# this is using equal rate transitions between all states (but you can change that with model = )
# this will print out which numbers map to which state
charP = char2PathsCategorical(phenvals, trees, useSpecies = allspecs,
                              plot = TRUE)

# or Generate paths using char2PathsCategorical (rooted, phytools version)
charP = char2PathsCategorical(phenvals, trees, useSpecies = allspecs,
                              use_rooted = TRUE, outgroup = "REFERENCE",
                              plot = TRUE)

# Correlate gene evolution with categorical trait using KW/Dunn (default)
cors = correlateWithCategoricalPhenotype(RERmat, charP)

# or Correlate gene evolution with categorical trait using ANOVA/Tukey
cors = correlateWithCategoricalPhenotype(RERmat, charP, method = "aov")

# Extract Results
# cors object (output of correlateWithCategoricalPhenotype) has two lists
# the first list is a data frame of the genes with Rho, P, and p.adj values
# the second list is a list of data frames, one for each pairwise analysis
# each data frame in the second list has two columns (Rho and P)
  # Rho is not the actual name of the statistic, but it is used to be consistent with the getStat function
allresults = cors[[1]]
pairwise_tables = cors[[2]]

# Enrichment Analysis
  # same as before, the res object is either "allresults" or cors[[1]]
  # or it is any of the tables in "pairwise_tables" i.e. pairwise_tables[[1]]

# read in the annotations
annots = readRDS("allannots.rds")

# run fastwilcoxGMTall
enrichments = fastwilcoxGMTall(getStat(allresults), annots, outputGeneVals=T)

# the third table is for the pairwise comparison between nocturnal and diurnal species
di_noc_enrichments = fastwilcoxGMTall(getStat(pairwise_tables[[3]]),
                                      annots, outputGeneVals=T)






