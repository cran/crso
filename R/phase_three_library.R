
# Exported:
### 1. makePhaseThreeImList(D,Q,rm.ordered,til.ee,pool.sizes,max.stored,max.nrs.borrow)

# Not exported
### 2. makeFullD1.IM(im.d1.core,borrow.pool,K)
### 3. makeFullD2.IM(im.d2.core,borrow.pool,K)
### 4. makeFullD3.IM(im.d3.core,borrow.pool,K)

### 5 updateTopIM.K3(D,Q.mat,rm,pool.sizes,max.nrs,til,max.stored)
### 6. updateTopIM(D,Q.mat,rm,pool.sizes,max.nrs,til,max.stored,K)



#' @title Make phase 3 im list from phase 2 im list
#'
#' @param D binary matrix of events by samples
#' @param Q penalty matrix of events by samples
#' @param rm.ordered matrix of rules ordered by phase one
#' @param til.ee list of rule set matrices (im list) from phase two
#' @param pool.sizes pool sizes for phase two
#' @param max.stored max number of rule sets saved
#' @param max.nrs.borrow max number of new rule sets per k, default is 10^5
#' @param shouldPrint Print progress updates? Default is TRUE
#' @examples
#' library(crso)
#' data(skcm)
#' list2env(skcm.list,envir=globalenv())
#' Q <- log10(P)
#' rm.full <- buildRuleLibrary(D,rule.thresh = 0.05) # Rule library matrix, dimension: 60 x 71
#' til.p2 <- makePhaseTwoImList(D,Q,rm.full,k.max = 3,pool.sizes=c(60,10,10),
#'           max.stored=100,shouldPrint = FALSE)
#' til.p3 <- makePhaseThreeImList(D,Q,rm.ordered = rm.full,til.ee = til.p2, pool.sizes=c(60,20,20),
#'          max.stored=100,max.nrs.borrow=100,shouldPrint = TRUE)
#' @export
#' @return phase 3 top im list
makePhaseThreeImList <- function(D,Q,rm.ordered,til.ee,pool.sizes,max.stored,max.nrs.borrow,shouldPrint){
  if(missing(shouldPrint)) shouldPrint <- TRUE
  if(missing(max.nrs.borrow)) max.nrs.borrow <- 10^5
  if(shouldPrint) print("Starting Phase 3: neighbor expansion")
  fm.ordered <- getFamMat(rm.ordered)
  ### Step 3: make updated top im lists
  grand.beg <- Sys.time()
  til.final <- til.ee ### til = top im list
  for(K in 3:length(til.final)){
    beg <- Sys.time()
    til.final[[K]] <- updateTopIM(D,Q,rm.ordered,pool.sizes,max.nrs.borrow,til.final,max.stored,K)
    if(shouldPrint) print(paste0("Updated top im for K = ",K,", time = ",signif(difftime(Sys.time(),beg,units="min"),4)," min"))
  }
  if(shouldPrint) print(paste0("Total Phase 3 Time: " ,signif(difftime(Sys.time(),grand.beg,units="min"),4)," min"))
  return(til.final)
}

### 2. This function takes in a matrix of unique cores of length K-1, and make IM of adding each borrow rule to every core
makeFullD1.IM <- function(im.d1.core,borrow.pool,K){
  full.d1.im <- matrix(NA,nrow=0,ncol=K)
  for(j in 1:nrow(im.d1.core)){
    cur.core <- im.d1.core[j,]
    core.im <- matrix(NA,nrow=length(borrow.pool),ncol=K)
    for(z in 1:nrow(core.im)) core.im[z,] <- c(cur.core,borrow.pool[z])
    full.d1.im <- rbind(full.d1.im,core.im)
  }
  return(full.d1.im)
}

### 3. Takes a core matrix of K-2 cores, and a borrow pool and makes all distance two rule sets

#' @importFrom utils combn
makeFullD2.IM <- function(im.d2.core,borrow.pool,K){
  full.d2.im <- matrix(NA,nrow=0,ncol=K)
  partial.borrow.im <- t(combn(borrow.pool,2))
  for(j in 1:nrow(im.d2.core)){
    cur.core <- im.d2.core[j,]
    core.im <- matrix(NA,nrow=nrow(partial.borrow.im),ncol=K)
    for(z in 1:nrow(core.im)) core.im[z,] <- c(cur.core,partial.borrow.im[z,])
    full.d2.im <- rbind(full.d2.im,core.im)
  }
  return(full.d2.im)
}

### 4. Takes a core matrix of K-3 cores, and a borrow pool and makes all distance 3 rule sets

#' @importFrom utils combn
makeFullD3.IM <- function(im.d3.core,borrow.pool,K){
  full.d3.im <- matrix(NA,nrow=0,ncol=K)
  partial.borrow.im <- t(combn(borrow.pool,3))
  for(j in 1:nrow(im.d3.core)){
    cur.core <- im.d3.core[j,]
    core.im <- matrix(NA,nrow=nrow(partial.borrow.im),ncol=K)
    for(z in 1:nrow(core.im)) core.im[z,] <- c(cur.core,partial.borrow.im[z,])
    full.d3.im <- rbind(full.d3.im,core.im)
  }
  return(full.d3.im)
}

### 5. Special handling of K = 3

#' @importFrom utils combn
updateTopIM.K3 <- function(D,Q.mat,rm,pool.sizes,max.nrs,til,max.stored){
  K <- 3
  top.im <- til[[K]]
  prev.im <- til[[K-1]]
  if(nrow(rm) <= pool.sizes[K]) return(top.im) ### if all rules have been exhaustively analyzed, should alwasy be = not less
  rm.cm <- makeRSCoverageMat(D,rm)
  first.new.rule <- pool.sizes[K] + 1
  cur.top.rs <- top.im[1,] ### Length K
  prev.top.rs <- prev.im[1,] ### Length K-1

  ### Special handling, only do this for K = 3:
  all.best.rs <- union(cur.top.rs,prev.top.rs)
  full.borrow.pool <- c(first.new.rule:nrow(rm)) ### always can do full for dist1

  ### Get dist 1 top im:
  im.d1.core <- t(combn(all.best.rs,K-1))   ### Let's get all k-1 cores
  full.d1.im <- makeFullD1.IM(im.d1.core,full.borrow.pool,K)
  top.im.d1 <- getTopIm(D,rm,rm.cm,full.d1.im,Q.mat,max.stored,should.check = TRUE)

  ### Get dist 2 top im:
  im.d2.core <- t(combn(all.best.rs,K-2))   ### Let's get all k-1 cores
  max.nrs <- max.nrs/nrow(im.d2.core) ### Adjust max nrs for number of cores
  ### Get borrow pool based on max nrs
  bp.length <- min(getMaxNforK(2,max.nrs),length(full.borrow.pool))
  borrow.pool <- full.borrow.pool[1:bp.length]
  full.d2.im <- makeFullD2.IM(im.d2.core,borrow.pool,K)
  top.im.d2 <- getTopIm(D,rm,rm.cm,full.d2.im,Q.mat,max.stored,should.check = TRUE)

  new.top.im <- rbind(top.im,top.im.d1,top.im.d2)
  new.top.im <- getTopIm(D,rm,rm.cm,new.top.im,Q.mat,max.stored,should.check = FALSE)
  return(new.top.im)
}


### 6. For all K, special handling of K = 3
updateTopIM <- function(D,Q.mat,rm,pool.sizes,max.nrs,til,max.stored,K){
  if(K==3) return(updateTopIM.K3(D,Q.mat,rm,pool.sizes,max.nrs,til,max.stored))
  top.im <- til[[K]]
  prev.im <- til[[K-1]]
  if(nrow(rm) <= pool.sizes[K]) return(top.im) ### if all rules have been exhaustively analyzed, should alwasy be = not less
  rm.cm <- makeRSCoverageMat(D,rm)
  first.new.rule <- pool.sizes[K] + 1
  cur.top.rs <- top.im[1,] ### Length K
  prev.top.rs <- prev.im[1,] ### Length K-1
  full.borrow.pool <- c(first.new.rule:nrow(rm)) ### always can do full for dist1

  ### Make D1 core:
  im.d1.core <- t(combn(cur.top.rs,K-1))
  im.d1.core <- rbind(im.d1.core,prev.top.rs)
  im.d1.core <- im.d1.core[!duplicated(im.d1.core),]
  ### Get dist 1 top im:
  full.d1.im <- makeFullD1.IM(im.d1.core,full.borrow.pool,K)
  top.im.d1 <- getTopIm(D,rm,rm.cm,full.d1.im,Q.mat,max.stored,should.check = TRUE)

  ### Get dist 2 top im:
  im.d2.core <- t(combn(cur.top.rs,K-2))   ### Let's get all k-1 cores
  im.d2.core <- rbind(im.d2.core,t(combn(prev.top.rs,K-2)))
  im.d2.core <- im.d2.core[!duplicated(im.d2.core),]
  max.nrs.d2 <- max.nrs/nrow(im.d2.core) ### Adjust max nrs for number of cores
  ### Get borrow pool based on max nrs
  bp.length <- min(getMaxNforK(2,max.nrs.d2),length(full.borrow.pool))
  borrow.pool <- full.borrow.pool[1:bp.length]
  #borrow.pool <- full.borrow.pool[1:getMaxNforK(2,max.nrs.d2)]
  full.d2.im <- makeFullD2.IM(im.d2.core,borrow.pool,K)
  top.im.d2 <- getTopIm(D,rm,rm.cm,full.d2.im,Q.mat,max.stored,should.check = TRUE)


  ### Get dist 3 top im:
  im.d3.core <- t(combn(cur.top.rs,K-3))   ### Let's get all k-1 cores
  im.d3.core <- rbind(im.d3.core,t(combn(prev.top.rs,K-3)))
  im.d3.core <- im.d3.core[!duplicated(im.d3.core),]
  if(K==4) im.d3.core <- as.matrix(im.d3.core) ### because duplicated makes it vector
  max.nrs.d3 <- max.nrs/nrow(im.d3.core) ### Adjust max nrs for number of cores
  ### Get borrow pool based on max nrs
  bp.length <- min(getMaxNforK(3,max.nrs.d3),length(full.borrow.pool))
  borrow.pool <- full.borrow.pool[1:bp.length]
  #borrow.pool <- full.borrow.pool[1:getMaxNforK(2,max.nrs.d3)]
  full.d3.im <- makeFullD3.IM(im.d3.core,borrow.pool,K)
  top.im.d3 <- getTopIm(D,rm,rm.cm,full.d3.im,Q.mat,max.stored,should.check = TRUE)


  new.top.im <- rbind(top.im,top.im.d1,top.im.d2,top.im.d3)
  new.top.im <- getTopIm(D,rm,rm.cm,new.top.im,Q.mat,max.stored,should.check = FALSE)
  return(new.top.im)
}


