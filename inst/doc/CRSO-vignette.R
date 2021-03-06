## ---- include = FALSE----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup---------------------------------------------------------------
library(crso)
### Load example dataset consisting of TCGA melanoma (SKCM) patients.
data(skcm)
list2env(skcm.list,envir=globalenv())
names(skcm.list) ### load D, P and cnv.dictionary
Q <- log10(P) ### Q is the penalty matrix derived from P

## ------------------------------------------------------------------------
rule.thresh <- .05 # Minimum percentage of samples covered by each rule in the rule library. (Default = .03)
spr <- 4 # Phase 1 random sets per rule (default = 40, recommend at least 40) 
trn <- 40 # Phase 1 stop criteria (default = 16, recommend at most 24)
k.max <- 6 # Maximum RS size (default is 12)
max.stored <- 1000 # Max RS stored per K in phases 2 and 3 (no default, recommend at least 10^4)
max.nrs.ee <- 1000 # Phase 2 max rs per K (default is 10^5, recommend at least 10^5)
max.compute <- 10*max.nrs.ee # Max raw im size for phase 2 (default is 10^6, recommend at least 10*max.nrs.ee)
max.nrs.borrow <- 1000 # Phase 3 max rs per K (default is 10^5, recommend at least 10^5)
filter.thresh <- .03 # minimum assignment threshold per rule set (default is .03)
### Generalized core parameters:
num.gc.iter <- 10 # Number GC iterations (default is 100)
num.gc.eval <- 100 # Rulesets evaluated per K per GC iter (default is 1000)

## ----phase 1-------------------------------------------------------------
set.seed(100)
rm.full <- buildRuleLibrary(D,rule.thresh) # build rule library
rm.ordered <- makePhaseOneOrderedRM(D,rm.full,spr,Q,trn,shouldPrint = TRUE) # run phase 1

## ----phase 2-------------------------------------------------------------
pool.sizes <- getPoolSizes(rm.ordered,k.max,max.nrs.ee,max.compute) 
### The pool size for each K is the number of rules considered for exhaustive evaluation in phase 2.
til.p2 <- makePhaseTwoImList(D,Q,rm.ordered,k.max,pool.sizes,max.stored,shouldPrint = TRUE) # Run phase 2
### TIL stands for top index list. The output of phase two is a list of top index matrices for each k.  Each index matrix contains the rule sets ordered by performance. For example the best performing rule set of size 3 will be the first row of the K.3 index matrix. For K=1 the index matrix is actually a vector.

## ----phase 3-------------------------------------------------------------
til.p3 <- makePhaseThreeImList(D,Q,rm.ordered,til.p2,pool.sizes,max.stored,max.nrs.borrow,shouldPrint = TRUE)
### Make TIL for phase 3 by updating phase two til to consider neighbor rule sets.

til.filtered <- makeFilteredImList(D,Q,rm.ordered,til.p3,filter.thresh)
### Filter the phase 3 results to only include rule sets for which every rule is assigned to a minimum percentage of samples, default is 3%
tpl.filtered <- evaluateListOfIMs(D,Q,rm.ordered,til.filtered) 
### Get top performance list (TPL), which contains the objective function score of all of the rule sets in til.filtered

## ----get core------------------------------------------------------------
best.rs.list <- getBestRsList(rm.ordered,tpl.filtered,til.filtered)
### This is a list of the best rule sets for all K

core.K <- getCoreK(D,rm.ordered,tpl.filtered,til.filtered) # Determine core K
core.ruleset <- getCoreRS(D,rm.ordered,tpl.filtered,til.filtered) # Extract core rule set
print(core.ruleset)

## ----gen core------------------------------------------------------------
list.subset.cores <- makeSubCoreList(D,Q,rm.ordered,til.filtered,num.gc.iter,num.gc.eval)
### list.subset.cores is a list of core rule set derived from subsampling iterations

gcr.df <- getGCRs(list.subset.cores) # Generalized core rules
print(gcr.df)
gcd.df <- getGCDs(list.subset.cores) # Generalized core duos
print(gcd.df)
gce.df <- getGCEs(list.subset.cores) # Generalized core events
print(gce.df)

