% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RERfuncs.R
\name{foreground2Tree}
\alias{foreground2Tree}
\title{Creates a binary trait tree from a set of foreground species.}
\usage{
foreground2Tree(
  foreground,
  treesObj,
  plotTree = T,
  clade = c("ancestral", "terminal", "all"),
  weighted = F,
  transition = "unidirectional",
  useSpecies = NULL
)
}
\arguments{
\item{treesObj}{A treesObj created by \code{\link{readTrees}}}

\item{plotTree}{Plot a tree representation of the result}

\item{clade}{A character string indicating which branches within the clade
containing the foreground species should be set to foreground. Must be one
of the strings "ancestral", "terminal", "all".}

\item{weighted}{if set to TRUE weights foreground edges belonging to the same clade such that their branch lengths sum up to 1 (only done for clade options "all" and "terminal").}

\item{transition}{A character string indicating whether transitions between background and foreground branches
are "bidirectional" or "unidirectional" (no foreground to background transitions, the default)}

\item{useSpecies}{Give only a subset of the species to use for ancestral state reconstruction
(e.g., only those species for which the trait can be reliably determined).}

\item{foreground.}{A character vector containing the foreground species}

\item{collapse2anc}{Put all the weight on the ancestral branch when the trait appears on a whole clade
(redundant to "clade", kept for backwards compatibility)}

\item{wholeClade}{Whether to implement the weighted edge option across
all members of a foreground clade (redundant to "clade", kept for backwards compatibility)}
}
\value{
A tree with edge.lengths representing phenotypic states
}
\description{
Creates a binary trait tree from a set of foreground species.
}
