
### This file will contain functions for making GCRs, GCDs and GCEs

### 1. oneSubCoreIteration(D,Q,rm,tpl,til,subset.size,num.evaluated)
### 2. makeSubCoreList(D,Q,rm,tpl,til,num.subsets,num.evaluated) - Exported
### 3. makeConfLevel(conf)
### 4. getGCRs(list.subsets.cores) - Exported
### 5. getDuosOneRule(rule)
### 6. getDuos(rules) # Get unique duos from multiple rules
### 7. getGCDs(list.subsets.cores) - Exported
### 8. getRuleEvents(rule)
### 9. getGCEs(list.subsets.cores) - Exported


# subset.size <- 0.8
#
# idc <- sample.int(ncol(D),floor(ncol(D)*subset.size),replace = FALSE)
# D.sub <- D[,idc]
# Q.sub <- Q[,idc]
# rm.cm.sub <- makeRSCoverageMat(D.sub,rm)
# sub.rs.idc.list <- getSubBestRsIdcList(D.sub,Q.sub,rm.imp,tpl,til,num.evaluated)

### 1.

#' @importFrom stats runif
oneSubCoreIteration <- function(D,Q,rm,til,subset.size,num.evaluated){
  # if(missing(subset.size)) subset.size <- 0.8
  # if(missing(num.evaluated)) num.evaluated <- 100

  idc <- sample.int(ncol(D),floor(ncol(D)*subset.size),replace = FALSE)
  D.sub <- D[,idc]
  Q.sub <- Q[,idc]
  rm.cm.sub <- makeRSCoverageMat(D.sub,rm)

  sub.rs.idc.list <- vector("list",length=length(til))
  names(sub.rs.idc.list) <- names(til)
  K <- 1
  sub.rs.idc.list[[K]] <- which.max(getSingleRuleWs(D.sub,rm,Q.sub))
  mini.beg <- Sys.time()
  for(K in 2:length(til)){
    im <- til[[K]]
    if(nrow(t(im))!=1) im <- im[1:min(nrow(im),num.evaluated),]

    if(nrow(t(im))==1) im <- as.matrix(t(im))
    perfs <- evaluateIM(D.sub,rm,rm.cm.sub,im,Q.sub)
    sub.rs.idc.list[[K]] <- im[which.max(perfs),]
  }
  core.cov.thresh <- runif(1,min=85,max=98)
  core.perf.thresh <- runif(1,min=80,max=95)

  ### Get best perfs and best covs
  best.sub.perfs <- rep(NA,length=length(sub.rs.idc.list))
  best.sub.covs <- best.sub.perfs
  for(K in 1:length(sub.rs.idc.list)){
    rs.idc <- sub.rs.idc.list[[K]]
    best.sub.perfs[K] <- getWofRS(D.sub,rm,rm.cm.sub,rs.idc,Q.sub)
    if(K==1) best.sub.covs[K] <- mean(rm.cm.sub[rs.idc,])
    if(K>1) best.sub.covs[K] <- mean(colSums(rm.cm.sub[rs.idc,])>0)
  }

  best.perf.percents <- 100*best.sub.perfs/max(best.sub.perfs)
  best.covs.percents <- 100*best.sub.covs/max(best.sub.covs)
  ### Determine core K
  core.idc <- which(best.covs.percents>=core.cov.thresh)
  core.idc <- intersect(core.idc,which(best.perf.percents>=core.perf.thresh))
  core.K <- min(core.idc)
  if(length(core.idc)==0) {
    #print("no K satisfies core")
    core.K <- length(best.sub.perfs)
  }
  sub.core.rules <- getRulesAsStrings(rm[sub.rs.idc.list[[core.K]],])
  return(sub.core.rules)
}

### 2.
#' @title Get list of core rules from random subsets of samples
#'
#' @param D input matrix D
#' @param Q input matrix Q
#' @param rm binary rule matrix
#' @param til list of top rule set index matrices
#' @param num.subsets number of subset iterations, default is 100
#' @param num.evaluated number of top rs considered per k per iteration, default is 1000
#' @param shouldPrint Print progress updates? Default is TRUE
#' @examples
#' library(crso)
#' data(skcm)
#' list2env(skcm.list,envir=globalenv())
#' Q <- log10(P)
#' rm.full <- buildRuleLibrary(D,rule.thresh = 0.05) # Rule library matrix, dimension: 60 x 71
#' til.p2 <- makePhaseTwoImList(D,Q,rm.full,k.max = 3,
#'           pool.sizes=c(60,20,20),max.stored=100,shouldPrint = FALSE)
#' subcore.list <- makeSubCoreList(D,Q,rm.full,til.p2,num.subsets=3,num.evaluated=50)
#' @export
makeSubCoreList <- function(D,Q,rm,til,num.subsets,num.evaluated,shouldPrint){
  if(missing(num.subsets)) num.subsets <- 100
  if(missing(num.evaluated)) num.evaluated <- 1000
  if(missing(shouldPrint)) shouldPrint <- TRUE
  subset.size <- 0.8
  k.max <- length(til)
  list.subsets.cores <- vector("list",length=num.subsets)

  beg <- Sys.time()
  for(j in 1:num.subsets){
    list.subsets.cores[[j]] <- oneSubCoreIteration(D,Q,rm,til,subset.size,num.evaluated)
    if(shouldPrint) print(paste0("Subset Core Iteration = ",j))
  }
  if(shouldPrint) print(Sys.time()-beg)
  return(list.subsets.cores)
}


### 3. Decide confidence level thresholds in here
makeConfLevel <- function(conf){
  med.thresh <- 40
  high.thresh <- 80
  conf.level <- rep("Low",length(conf))
  conf.level[which(conf>=med.thresh)] <- "Medium"#paste0("Intermediate: [",med.thresh,", ",high.thresh,")")
  conf.level[which(conf>=high.thresh)] <- "High"#paste0("High: >= ",high.thresh)
  return(conf.level)
}

### 4
#' @title Get Generalized Core Rules
#'
#' @param list.subset.cores list of subset cores
#' @examples
#' list.subset.cores <- list(c("A.B.C","D.E","A.D"),c("A.C","B.C.D","D.E"),
#' c("A.B.C","D.E"),c("A.B.C","D.E","B.C.D"))
#' getGCRs(list.subset.cores) # Confidence column should be 100, 75, 50, 25, 25
#' @export
getGCRs <- function(list.subset.cores){
  conf <- 100*sort(table(unlist(list.subset.cores)))/length(list.subset.cores) ### gen core confidence level
  conf.level <- makeConfLevel(conf)
  rules <- names(conf) ### gen core rules
  ### Make df.r (r for rules)
  df.r <- data.frame(GCR=rules,Confidence=as.numeric(conf),Confidence.Level=conf.level)
  df.r$GCR <- factor(df.r$GCR, levels = df.r$GCR[order(df.r$Confidence)])
  df.r <- df.r[nrow(df.r):1,]
  return(df.r)
}



### 5. Get duos
getDuosOneRule <- function(rule){
  all.duos <- c()
  events <- strsplit(rule,"\\.")[[1]]
  if(length(events)==2) all.duos <- c(all.duos, paste0(events[1],".",events[2]))
  if(length(events)>2){
    duos.mat <- t(combn(events,2))
    for(j in 1:nrow(duos.mat)) all.duos <- c(all.duos,paste0(duos.mat[j,],collapse = "."))
  }
  all.duos <- sort(all.duos)
  return(all.duos)
}

### 6. Get unique duos from multiple rules
getDuos <- function(rules){
  all.duos <- c()
  for(j in 1:length(rules)){
    rule <- rules[j]
    all.duos <- c(all.duos,getDuosOneRule(rule))
  }
  all.duos <- sort(unique(all.duos))
  return(all.duos)
}

### 7.
#' @title Get Generalized Core Duos
#'
#' @param list.subset.cores list of subset cores
#' @examples
#' list.subset.cores <- list(c("A.B.C","D.E","A.D"),c("A.C","B.C.D","D.E"),
#' c("A.B.C","D.E"),c("A.B.C","D.E","B.C.D"))
#' getGCDs(list.subset.cores) # Confidence column should be 100, 100, 100, 75, 50, 25, 25
#' @export
getGCDs <- function(list.subset.cores){
  list.subset.duos <- list.subset.cores
  for(j in 1:length(list.subset.duos)) list.subset.duos[[j]] <- getDuos(list.subset.cores[[j]])
  temp <- getGCRs(list.subset.duos)
  colnames(temp)[1] <- "GCD"
  return(temp)
}

### 8.
getRuleEvents <- function(rule) strsplit(rule,"\\.")[[1]]


### 9.
#' @title Get Generalized Core Events
#'
#' @param list.subset.cores list of subset cores
#' @examples
#' list.subset.cores <- list(c("A.B.C","D.E","A.D"),
#' c("A.C","B.C.D","D.E"),c("A.B.C","D.E"),c("A.B.C","D.E","B.C.D"))
#' getGCEs(list.subset.cores) # Confidence column should be 100, 100, 100, 100, 100
#' @export
getGCEs <- function(list.subset.cores){
  list.subset.events <- list.subset.cores
  for(j in 1:length(list.subset.cores)) list.subset.events[[j]] <- unique(as.character(unlist(sapply(list.subset.cores[[j]],getRuleEvents))))
  temp <- getGCRs(list.subset.events)
  colnames(temp)[1] <- "GCE"
  return(temp)
}
