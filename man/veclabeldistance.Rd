% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unlabelledFuncs.R
\name{veclabeldistance}
\alias{veclabeldistance}
\title{Get D2 distance between two vectors of labels}
\usage{
veclabeldistance(lab1, lab2, samesize = TRUE, eps = 0, nottips = FALSE)
}
\arguments{
\item{lab1}{A vector of tree labels (as created by treelabels)}

\item{lab2}{A vector of tree labels (as created by treelabels)}

\item{samesize}{Logical: are the trees the same size?}

\item{eps:}{if the trees are not the same size (samesize=F), 
the components in the distance will be divided by the number of tips. 
Then eps*(difference in tip number) is added to preserve the metric property.}

\item{nottips:}{exclude the number of 1s from the distance (default F)}
}
\value{
D2 distance: the L2 norm comparing the vectors of labels
}
\description{
Get D2 distance between two vectors of labels
}
\examples{
veclabeldistance(treelabels(rtree(10)), treelabels(rtree(12))
}

