% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/makePoolSizes.R
\name{getPoolSizes}
\alias{getPoolSizes}
\title{Get pool sizes for phase 2}
\usage{
getPoolSizes(rm.ordered, k.max, max.nrs.ee, max.compute)
}
\arguments{
\item{rm.ordered}{binary rule matrix ordered from phase 1}

\item{k.max}{maximum rule set size}

\item{max.nrs.ee}{max number of rule sets per k}

\item{max.compute}{maximum raw rule sets considered per k}
}
\description{
Get pool sizes for phase 2
}
\examples{
library(crso)
data(skcm)
list2env(skcm.list,envir=globalenv())
rm.full <- buildRuleLibrary(D,rule.thresh = 0.05) # Rule library matrix, dimension: 60 x 71
rm.ordered <- rm.full # Skip phase one in this example
getPoolSizes(rm.ordered,k.max = 7,max.nrs.ee = 10000)
# [1] 60  60  40  23  18  16  15
}
