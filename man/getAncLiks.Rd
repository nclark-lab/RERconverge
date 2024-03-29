% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RERfuncs.R
\name{getAncLiks}
\alias{getAncLiks}
\title{Returns ancestral likelihoods at each node. Based on code from ace in ape and fitMk/rerootingMethod in phytools}
\usage{
getAncLiks(tree, tipvals, Q = NULL, rate_model = "ER", root_prior = "auto")
}
\arguments{
\item{tree}{Object of class phylo that has been pruned to only contain the species in tipvals}

\item{tipvals}{The phenotype data for the extant species in the tree in the same order as tree$tip.label and mapped to integers}

\item{Q}{A transition matrix, if NULL the transition matrix is fit with fit_mk from castor package}

\item{rate_model}{The rate model to use for fitting the transition matrix if one is not provided}

\item{root_prior}{The root prior used when fitting the transition matrix if one is not provided}
}
\value{
The ancestral likelihoods at each node in the tree
}
\description{
Returns ancestral likelihoods at each node. Based on code from ace in ape and fitMk/rerootingMethod in phytools
}
