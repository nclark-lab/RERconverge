% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plottingFuncs.R
\name{plotTreeCategorical}
\alias{plotTreeCategorical}
\title{Plot a phenotype tree generated by \code{\link{char2TreeCategorical}} with the branches colored by category}
\usage{
plotTreeCategorical(
  tree,
  category_names = NULL,
  master = NULL,
  node_states = NULL
)
}
\arguments{
\item{tree}{The phenotype tree returned by \code{\link{char2TreeCategorical}}}

\item{category_names}{The names of the categories in order of the corresponding numerical labels}

\item{master}{If provided, the tree will be plotted using the branch lengths of that tree, coloring branches by the phenotype in the phenotype tree. The master tree and phenotype tree should have the same topology and the same species. If not provided, use.edge.length is set to FALSE.}

\item{node_states}{A vector of the phenotype state at each node in the tree (in order of the nodes). If provided, vertical lines representing nodes in the tree are colored by state. It must be numerical and the integers must have the same mapping to the category names as in the phenotype tree.}
}
\value{
A plot of a phylogenetic tree with branches colored by phenotype category.
}
\description{
Plot a phenotype tree generated by \code{\link{char2TreeCategorical}} with the branches colored by category
}
