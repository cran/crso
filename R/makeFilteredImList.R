

### 1. getMSA.oneRS(rm,rm.cm,rs.idc,Q)
### 2. getMSAofIM(rm,rm.cm,im,Q)
### 3. getMSAofImList(rm,rm.cm,til,Q)
### 4. makeFilteredImList(D,Q,rm,til,filter.thresh) -- exported


### Get Min Samples Assigned of one rule
getMSA.oneRS <- function(rm,rm.cm,rs.idc,Q){
  assign.vec <- makeOptimalAssignment(rm,rm.cm,rs.idc,Q)
  tab <- table(assign.vec[which(assign.vec!=0)])
  if(length(tab)<length(rs.idc)) return(0) ### Return 0 if not all rules are assigned
  return(min(tab))
}

#' @importFrom foreach foreach %dopar%
getMSAofIM <- function(rm,rm.cm,im,Q){

  im.cur <- NULL

  n.splits <- foreach::getDoParWorkers()
  im.list <- splitMatIntoList(im,n.splits)
  ac.vec <- foreach(im.cur = im.list,.combine = c,.export=ls(envir=globalenv())) %dopar% {
    ac.vec <- apply(im.cur,1,function(x)getMSA.oneRS(rm,rm.cm,x,Q))
  }
  return(ac.vec)
}

getMSAofImList <- function(D,Q,rm,til){
  rm.cm <- makeRSCoverageMat(D,rm)
  msa.list <- til
  msa.list[[1]] <- rowSums(rm.cm)
  for(K in 2:length(til)){
    im <- til[[K]]
    msa.list[[K]] <- getMSAofIM(rm,rm.cm,im,Q)
    #print(paste0("k = ",K))
  }
  return(msa.list)
}

#' @title Make filtered im list from phase 3 im list
#'
#' @param D binary matrix of events by samples
#' @param Q penalty matrix of events by samples
#' @param rm matrix of rules ordered by phase one
#' @param til im list from phase 3
#' @param filter.thresh minimum percentage of samples assigned to each rule in rs
#' @examples
#' library(crso)
#' data(skcm)
#' list2env(skcm.list,envir=globalenv())
#' Q <- log10(P)
#' rm.full <- buildRuleLibrary(D,rule.thresh = 0.05) # Rule library matrix, dimension: 60 x 71
#' til.p2 <- makePhaseTwoImList(D,Q,rm.full,k.max = 3,
#'           pool.sizes=c(60,20,20),max.stored=100,shouldPrint = FALSE)
#' filtered.im.list <- makeFilteredImList(D,Q,rm.full,til.p2,filter.thresh = 0.05)
#' @export
#' @return filtered top im list
makeFilteredImList <- function(D,Q,rm,til,filter.thresh){
  #if(missing(cov.thresh)) cov.thresh <- .02
  if(missing(filter.thresh)) filter.thresh <- .03
  min.samp.thresh <- 5
  samp.thresh <- max(ceiling(filter.thresh*ncol(D)),min.samp.thresh)
  msa.list <- getMSAofImList(D,Q,rm,til)

  new.im.list <- vector("list",length=length(til))
  names(new.im.list) <- paste0("K.",1:length(til))
  new.im.list[[1]] <- til[[1]][1]

  til.filtered <- vector("list",length=length(til))
  names(til.filtered) <- paste0("K.",1:length(til))
  #tpl.filtered <- til.filtered
  idc <- which(msa.list[[1]]>=samp.thresh)
  til.filtered[[1]] <- til[[1]][idc]
  #tpl.filtered[[1]] <- tpl[[1]][idc]
  for(K in 2:length(til.filtered)){
    idc <- which(msa.list[[K]]>=samp.thresh)
    if(length(idc)>0){
      til.filtered[[K]] <- til[[K]][idc,]
      #tpl.filtered[[K]] <- tpl[[K]][idc]
    }
  }
  ks.filtered <- which(unlist(lapply(til.filtered,length))>0)
  til.filtered <- til.filtered[ks.filtered]
  #tpl.filtered <- tpl.filtered[ks.filtered]
  return(til.filtered)
}




#
#
# makeFilteredTIL <- function(rm.ordered,cov.thresh,min.samp.thresh){
#   if(missing(cov.thresh)) cov.thresh <- .02
#   if(missing(min.samp.thresh)) min.samp.thresh <- 6
#   #list2env(p2.results.list,envir=globalenv())
#   #rm.imp <- rm.ordered ### importance rm
#   rownames(rm.imp) <- paste0("r",1:nrow(rm.imp))
#   rm.cm.imp <- makeRSCoverageMat(D,rm.imp)
#
#   min.samps.assigned <- max(ceiling(cov.thresh*ncol(D)),min.samp.thresh)
#   new.rs.idc.list <- vector("list",length=length(tpl.final))
#   names(new.rs.idc.list) <- paste0("K.",1:length(tpl.final))
#   new.rs.idc.list[[1]] <- til.final[[1]][1]
#
#   all.idcs <- rep(NA,length(til.final))
#   for(K in 2:length(all.idcs)){
#     im <- til.final[[K]]
#     go <- 1
#     j <- 1
#     while(go){
#       rs.idc <- im[j,]
#       assign.vec <- makeOptimalAssignment(rm.imp,rm.cm.imp,rs.idc,Q.mat)
#       assign.vec <- assign.vec[which(assign.vec!=0)]
#       if(min(table(assign.vec))>=min.samps.assigned) {
#         go <- 0
#         idc <- j
#       }
#       j <- j + 1
#       if(j%%500==0)print(j)
#       if(j>nrow(im)){
#         go <- 0
#         idc <- -1
#       }
#     }
#     all.idcs[K] <- idc
#     #print(all.idcs)
#     if(idc >= 1) new.rs.idc.list[[K]] <- im[idc,]
#     if(idc==-1) return(new.rs.idc.list)
#   }
#   return(new.rs.idc.list)
# }
# ################################################################################
#
#
#
#
#
# beg <- Sys.time()
# msa.list <- getMSAofImList(D,rm.ordered,til.p3,tpl.p3,Q)
# print(Sys.time()-beg)
#
#
#
# list.msa.lists <- vector("list",length=length(our.tissues))
# names(list.msa.lists) <- our.tissues
# beg <- Sys.time()
# for(tissue in our.tissues){
#   print(tissue)
#   load(paste0(p2.result.dir,tissue,"_p2_results.RData"))
#   list2env(p2.results.list,envir=globalenv())
#   rm.cm.ordered <- makeRSCoverageMat(D,rm.ordered)
#   list.msa.lists[[tissue]] <- getMSA.List(rm.ordered,rm.cm.ordered,til.final,tpl.final,Q=Q.mat)
#   print(Sys.time()-beg)
# }
# outfile <- "PAN.TISSUE.FAM/all_msa_lists.RData"
# save(list.msa.lists,file=outfile)
# makeFullResultList <- function(p2.results.list,filter.thresh,min.samp.thresh,msa.list){
#   list2env(p2.results.list,envir=globalenv())
#
# }
#
# ###################################################################################
# ### Main
# ###################################################################################
# tissue <- "SKCM"
# load(paste0(p2.result.dir,tissue,"_p2_results.RData"))
# msa.list <- list.msa.list[[tissue]]
#
#
# ###################################################################################
# assign.vec <- makeOptimalAssignment(rm.imp,rm.cm.imp,rs.idc,Q.mat)
# assign.vec <- assign.vec[which(assign.vec!=0)]
#
# makeFilteredTIL <- function(p2.results.list,cov.thresh,min.samp.thresh){
#   if(missing(cov.thresh)) cov.thresh <- .02
#   if(missing(min.samp.thresh)) min.samp.thresh <- 6
#   list2env(p2.results.list,envir=globalenv())
#   rm.imp <- rm.ordered ### importance rm
#   rownames(rm.imp) <- paste0("r",1:nrow(rm.imp))
#   rm.cm.imp <- makeRSCoverageMat(D,rm.imp)
#
#   min.samps.assigned <- max(ceiling(cov.thresh*ncol(D)),min.samp.thresh)
#   new.rs.idc.list <- vector("list",length=length(tpl.final))
#   names(new.rs.idc.list) <- paste0("K.",1:length(tpl.final))
#   new.rs.idc.list[[1]] <- til.final[[1]][1]
#
#   all.idcs <- rep(NA,length(til.final))
#   for(K in 2:length(all.idcs)){
#     im <- til.final[[K]]
#     go <- 1
#     j <- 1
#     while(go){
#       rs.idc <- im[j,]
#       assign.vec <- makeOptimalAssignment(rm.imp,rm.cm.imp,rs.idc,Q.mat)
#       assign.vec <- assign.vec[which(assign.vec!=0)]
#       if(min(table(assign.vec))>=min.samps.assigned) {
#         go <- 0
#         idc <- j
#       }
#       j <- j + 1
#       if(j%%500==0)print(j)
#       if(j>nrow(im)){
#         go <- 0
#         idc <- -1
#       }
#     }
#     all.idcs[K] <- idc
#     #print(all.idcs)
#     if(idc >= 1) new.rs.idc.list[[K]] <- im[idc,]
#     if(idc==-1) return(new.rs.idc.list)
#   }
#   return(new.rs.idc.list)
# }
# ################################################################################




