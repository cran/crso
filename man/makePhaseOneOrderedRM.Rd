% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phase_one_library.R
\name{makePhaseOneOrderedRM}
\alias{makePhaseOneOrderedRM}
\title{Order rules according to phase one importance ranking}
\usage{
makePhaseOneOrderedRM(D, rm.start, spr, Q, trn, n.splits, shouldPrint)
}
\arguments{
\item{D}{Binary matrix of N events and M samples}

\item{rm.start}{Starting binary rule matrix (i.e., rule library)}

\item{spr}{Random rule sets per rule in each phase one iteration. Default is 40.}

\item{Q}{Penalty matrix, negative log of passenger probability matrix.}

\item{trn}{Target rule number for stopping iterating. Default is 16.}

\item{n.splits}{number of splits for parallelization. Default is all available cpus.}

\item{shouldPrint}{Print progress updates? Default is TRUE}
}
\value{
binary rule matrix ordered by phase one importance ranking
}
\description{
Order rules according to phase one importance ranking
}
\examples{
data(skcm)
list2env(skcm.list,envir=globalenv())
Q <- log10(P)
rm.full <- buildRuleLibrary(D,rule.thresh = 0.06) # Rule library matrix, dimension: 36s x 71
rm.ordered <- makePhaseOneOrderedRM(D,rm.full,spr = 1,Q,trn = 34,shouldPrint = TRUE)
# note, for real applications, spr should be at least 40.
}
