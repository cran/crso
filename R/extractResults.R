
### 1. getBestRsCoverages(D,rm,tpl,til)
### 2. getBestRsIdcList(rm,tpl,til)
### 3. getBestRsList(rm,tpl,til) - Exported
### 4. getCoreK(D,rm,tpl,til,cov.thresh,perf.thresh) - Exported
### 5. getCoreRS(D,rm,tpl,til,cov.thresh,perf.thresh) - Exported


############################################################################################
# 1.
getBestRsCoverages <- function(D,rm,tpl,til){
  best.coverages <- rep(NA,length(tpl))
  rm.cm <- makeRSCoverageMat(D,rm)
  best.rs.idc.list <- getBestRsIdcList(rm,tpl,til)
  best.coverages[1] <- mean(rm.cm[best.rs.idc.list[[1]],])
  for(k in 2:length(tpl)){
    rs.idc <- best.rs.idc.list[[k]]
    best.coverages[k] <- length(which(colSums(rm.cm[rs.idc,])>0))/ncol(D)
  }
  return(best.coverages)
}



# 2. Get best rs.idc of each size
getBestRsIdcList <- function(rm,tpl,til){
  best.rs.list <- tpl ### Initialization
  for(j in 1:length(tpl)){
    cur.perfs <- tpl[[j]]
    cur.im <- til[[j]]
    if(length(cur.perfs)==1) rs.idc <- cur.im
    if(length(cur.perfs)>1){
      if(j>1)rs.idc <- cur.im[which.max(cur.perfs),]
      if(j==1) rs.idc <- cur.im[which.max(cur.perfs)]
    }
    best.rs.list[[j]] <- rs.idc #getRulesAsStrings(rm[rs.idc,])
  }
  return(best.rs.list)
}

# 3. ### Get best rule set of each size
#' @title Get list of best rule sets of size K for all K
#'
#' @param rm binary rule matrix
#' @param tpl list of top performances
#' @param til list of top rule set index matrices
#' @examples
#' library(crso)
#' data(skcm)
#' list2env(skcm.list,envir=globalenv())
#' Q <- log10(P)
#' rm.full <- buildRuleLibrary(D,rule.thresh = 0.05) # Rule library matrix, dimension: 60 x 71
#' til.p2 <- makePhaseTwoImList(D,Q,rm.full,k.max = 3,pool.sizes=c(60,20,20),
#'           max.stored=100,shouldPrint = FALSE)
#' tpl.p2 <- evaluateListOfIMs(D,Q,rm.full,til.p2)
#' best.rs.list <- getBestRsList(rm = rm.full,tpl = tpl.p2,til = til.p2)
#' @export
getBestRsList <- function(rm,tpl,til){
  best.rs.list <- tpl ### Initialization
  for(j in 1:length(tpl)){
    cur.perfs <- tpl[[j]]
    cur.im <- til[[j]]
    if(length(cur.perfs)==1) rs.idc <- cur.im
    if(length(cur.perfs)>1){
      if(j>1)rs.idc <- cur.im[which.max(cur.perfs),]
      if(j==1) rs.idc <- cur.im[which.max(cur.perfs)]
    }
    best.rs.list[[j]] <- getRulesAsStrings(rm[rs.idc,])
  }
  return(best.rs.list)
}

### 4
#' @title Determine core K from phase 3 tpl and til
#'
#' @param D input matrix D
#' @param rm binary rule matrix
#' @param tpl list of top performances
#' @param til list of top rule set index matrices
#' @param cov.thresh core coverage threshold, defaults is 95
#' @param perf.thresh core performance threshold, default is 90
#' @examples
#' library(crso)
#' data(skcm)
#' list2env(skcm.list,envir=globalenv())
#' Q <- log10(P)
#' rm.full <- buildRuleLibrary(D,rule.thresh = 0.05) # Rule library matrix, dimension: 60 x 71
#' til.p2 <- makePhaseTwoImList(D,Q,rm.full,k.max = 3,pool.sizes=c(60,20,20),
#'           max.stored=100,shouldPrint = FALSE)
#' tpl.p2 <- evaluateListOfIMs(D,Q,rm.full,til.p2)
#' core.K <- getCoreK(D,rm.full,tpl.p2,til.p2)
#' # core.K should be 3 almost always for this example, can run a few time to confirm
#' @export
getCoreK <- function(D,rm,tpl,til,cov.thresh,perf.thresh){
  if(missing(cov.thresh)) cov.thresh <- 95
  if(missing(perf.thresh)) perf.thresh <- 90
  best.perfs <- unlist(lapply(tpl,max))
  best.covs <- getBestRsCoverages(D,rm,tpl,til)
  best.perf.percents <- 100*best.perfs/max(best.perfs)
  best.covs.percents <- 100*best.covs/max(best.covs)
  ### Determine core K
  core.idc <- which(best.covs.percents>=cov.thresh)
  core.idc <- intersect(core.idc,which(best.perf.percents>=perf.thresh))
  core.K <- min(core.idc)
  if(length(core.idc)==0) {
    #print("no K satisfies core")
    core.K <- length(best.perfs)
  }
  return(core.K)
}

### 5
#' @title Get core rules from phase 3 tpl and til
#'
#' @param D input matrix D
#' @param rm binary rule matrix
#' @param tpl list of top performances
#' @param til list of top rule set index matrices
#' @param cov.thresh core coverage threshold, defaults is 95
#' @param perf.thresh core performance threshold, default is 90
#' @examples
#' library(crso)
#' data(skcm)
#' list2env(skcm.list,envir=globalenv())
#' Q <- log10(P)
#' rm.full <- buildRuleLibrary(D,rule.thresh = 0.05) # Rule library matrix, dimension: 60 x 71
#' til.p2 <- makePhaseTwoImList(D,Q,rm.full,k.max = 3,pool.sizes=c(60,20,20),
#'           max.stored=100,shouldPrint = FALSE)
#' tpl.p2 <- evaluateListOfIMs(D,Q,rm.full,til.p2)
#' core.rs <- getCoreRS(D,rm.full,tpl.p2,til.p2) # core.rs should be r1, r2, r3
#' @export
getCoreRS <- function(D,rm,tpl,til,cov.thresh,perf.thresh){
  if(missing(cov.thresh)) cov.thresh <- 95
  if(missing(perf.thresh)) perf.thresh <- 90
  best.rs.list <- getBestRsList(rm,tpl,til)
  best.perfs <- unlist(lapply(tpl,max))
  best.covs <- getBestRsCoverages(D,rm,tpl,til)
  core.K <- getCoreK(D,rm,tpl,til,cov.thresh,perf.thresh)
  core.rules <- best.rs.list[[core.K]]
  return(core.rules)
}






