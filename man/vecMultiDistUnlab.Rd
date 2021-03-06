% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unlabelledFuncs.R
\name{vecMultiDistUnlab}
\alias{vecMultiDistUnlab}
\title{All pairwise D2 distances for a list of trees}
\usage{
vecMultiDistUnlab(trees, listoflabels = NA, samesize = TRUE, eps = 0,
  nottips = FALSE)
}
\arguments{
\item{trees}{A list of trees}

\item{listoflabels}{Alternatively, a list of labels. If listoflabels is NA
(the default), this function first computes the labels for the set of trees.
If listoflabels is not NA, the labels are compared directly and the trees are ignored.}

\item{samesize}{Logical: whether the trees are the same size}

\item{eps:}{if the trees are not the same size (samesize=F), the components 
in the distance will be divided by the number of tips.
Then eps*(difference in tip number) is added to preserve the metric property.}

\item{nottips:}{exclude the number of 1s from the distance (default F)}
}
\value{
An object of class 'dist' containing tree-tree D2 distances (uses veclabeldistance)
}
\description{
All pairwise D2 distances for a list of trees
}
\examples{
dd=vecmultiDistUnlab(rmtree(10,12))
}

