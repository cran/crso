### This file contains the functions for evaluating the performance of a rule set.
### W refers to the objective function score.

### None of the functions in this file are exported.

### 1) makeOptimalAssignment(rm,rm.cm,rs.idc,Q.mat)
### 2) makePassengerD(D,rm,rs.idc,assign.vec)
### 3) getWofAssignment(D,rm,rs.idc,assign.vec,Q.mat,W.0)
### 4) getWofRS(D,rm,rm.cm,rs.idc,Q.mat,W.0)
### 5) getSingleRuleWs(D,rm,Q.mat)




# 1
makeOptimalAssignment <- function(rm,rm.cm,rs.idc,Q.mat){
  if(length(rs.idc)==0) return(rep(0,ncol(rm.cm)))
  if(length(rs.idc)==1) {
    assign.vec <- as.numeric(rm.cm[rs.idc,])
    assign.vec[which(assign.vec==1)] <- rs.idc
    return(assign.vec)
  }
  assign.vec <- rep(0,length=ncol(rm.cm))
  rs.cm <- rm.cm[rs.idc,]
  for(j in 1:ncol(rs.cm)){
    candidate.rule.idc <- rs.idc[which(rs.cm[,j]==1)]
    if(length(candidate.rule.idc)==1) assign.vec[j] <-  candidate.rule.idc
    if(length(candidate.rule.idc)>1){
      ### here we score each rule
      rule.scores <- rep(NA,length=length(candidate.rule.idc))
      for(y in 1:length(rule.scores)){
        rule.id <- candidate.rule.idc[y]
        rule.events <- which(rm[rule.id,]==1)
        rule.scores[y] <- sum(Q.mat[rule.events,j])
      }
      assign.vec[j] <- candidate.rule.idc[which.min(rule.scores)]
    }
  }
  return(assign.vec)
}

# 2.
makePassengerD <- function(D,rm,rs.idc,assign.vec){
  #rs.mat <- rm[rs.idc,]
  passenger.D <- D
  ### if rs in no rules
  if(length(rs.idc)==0)return(D)
  ### if rs is only one rule
  if(length(rs.idc)==1){
    rule.events.idc <- which(rm[rs.idc,]==1)
    passenger.D[rule.events.idc,which(assign.vec==rs.idc)] <- 0
    return(passenger.D)
  }
  for(j in 1:ncol(D)){
    if(assign.vec[j]>0){
      rule <- assign.vec[j]
      rule.events.idc <- which(rm[rule,]==1)
      passenger.D[rule.events.idc,j] <- 0
    }
  }
  return(passenger.D)
}


# 3.
getWofAssignment <- function(D,rm,rs.idc,assign.vec,Q.mat,W.0){
  if(missing(W.0)) W.0 <- sum(Q.mat)
  ### We pass along Q.0 to avoid
  D.pass <- makePassengerD(D,rm,rs.idc,assign.vec)
  Q.assigned <- Q.mat
  Q.assigned[which(D.pass==0)] <- 0
  W <- sum(Q.assigned) - W.0
  return(W)
}


# 4.
getWofRS <- function(D,rm,rm.cm,rs.idc,Q.mat,W.0){
  if(missing(W.0)) W.0 <- sum(Q.mat)
  rs.idc <- unique(rs.idc)
  assign.vec <- makeOptimalAssignment(rm,rm.cm,rs.idc,Q.mat)
  #assign.vec <- makeHierarchyAssignment(rm,rm.cm,rs.idc)
  W <- getWofAssignment(D,rm,rs.idc,assign.vec,Q.mat,W.0)
  return(W)
}

# 5.
getSingleRuleWs <- function(D,rm,Q.mat){
  W.0 <- sum(Q.mat)
  perfs.vec <- rep(NA,nrow(rm))
  names(perfs.vec) <- rownames(rm)
  rm.cm <- makeRSCoverageMat(D,rm)
  for(rule.id in 1:nrow(rm))perfs.vec[rule.id] <- getWofRS(D,rm,rm.cm,rs.idc=rule.id,Q.mat,W.0)
  return(perfs.vec)
}



