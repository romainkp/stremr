###################################
## SE of IPW estimator of MSM coef
###################################
getSEcoef <- function(ID, clusterID = NULL, nID, clusterNperID = clusterNperID, t.var, Yname, MSMdata, MSMpredict, MSMdesign, IPW_MSMestimator = TRUE, family="binomial", stabilize, use_weights){
  cum.IPAW <- "cum.IPAW"
  
  if(is.null(clusterID))clusterID <- ID
  
  if(!is.null(stabilize) & use_weights){ #getSEcoef is called from survMSM (i.e. the MSM is saturated) & using weights from IPW
      ##IC formula derived from appendix of "A Case Study of the Impact of Data-Adaptive Versus Model-Based Estimation of the Propensity Scores on Causal Inferences from Three Inverse Probability Weighting Estimators" and after working out the IPW normalizing constant under the assumption that the action mechanism is estimated consistently and a saturated MSM for hazard
      
    saturatedMSMterms <- colnames(MSMdesign)[grep("Periods",colnames(MSMdesign))]
    IC.data <- MSMdata[ ,c(unique(c(ID,clusterID)), t.var, Yname, cum.IPAW, MSMpredict,"cum.stab.P",saturatedMSMterms), with = FALSE]
    
    ## Calculate S_theta(t-1) = P(Y_theta(0)=0, ..., Y_theta(t)=0) because Y_theta(t) = I( T_theta in interval t-1) where T_theta is failure time of interest
    Stheta <- unique(MSMdata[, c(cum.IPAW,t.var, "rule.name", MSMpredict, saturatedMSMterms), with = FALSE][cum.IPAW>0,][,cum.IPAW:=NULL])
    setkeyv(Stheta,c("rule.name",t.var))
    Stheta[,Stheta.t := cumprod(1-get(MSMpredict)), by=rule.name]
    Stheta[,Stheta.tminus1 := shift(Stheta.t, n=1L, type = "lag", fill = 1), by=rule.name]
    
    ## very important to set sort to FALSE here so that ordering of rows remain unchanged:
    IC.data <- merge(IC.data, Stheta[,c(saturatedMSMterms,"Stheta.tminus1"), with=FALSE], by = saturatedMSMterms, all.x=TRUE, all.y=FALSE,sort=FALSE)
    IC.data[, (saturatedMSMterms) :=NULL ]
    
    IC.data[, "ModifiedEpsilon" := ifelse( get(Yname)==1,  1/get(MSMpredict), -1/(1-get(MSMpredict)) )/Stheta.tminus1 ]
    IC.data[is.infinite(ModifiedEpsilon), ModifiedEpsilon:=.Machine$double.xmax] # replace Inf values with largest possible double value so that when multiplied by 0 it gives 0 and not NaN 
    
    if (IPW_MSMestimator) {
      if(stabilize){ # undo stabilization
        sweep.tmp <- (IC.data[[cum.IPAW]]/IC.data[["cum.stab.P"]]) * IC.data[["ModifiedEpsilon"]]
      }else{ # nothing to undo, use unstabilized weights
        sweep.tmp <- IC.data[[cum.IPAW]] * IC.data[["ModifiedEpsilon"]]
      }
    }else{
      sweep.tmp <- IC.data[["ModifiedEpsilon"]]
    }
    
    D.O.beta <- sweep(t(MSMdesign), 2, sweep.tmp, "*") # multiply each column of t(MSMdesign) by sweep.tmp:
    xDT <- data.table(ID = IC.data[[clusterID]], t(D.O.beta)) ## sum of t and theta with no clustering and over t, theta, and cluster observations otherwise
    setkeyv(xDT, cols = "ID")
    xDT <- xDT[, lapply(.SD, sum), by = ID] ## take sum by ID, no longer drop ID column
    cols <- setdiff(names(xDT), "ID")
    xDT[,(cols) := nID/clusterNperID[,sum(clusterN)]*.SD, .SDcols = cols] # apply the clustering weights (n/N where n is iid sample size and N = sum_cluster of cluster size
    D.O.beta.alt <- as.matrix(xDT, rownames = "ID") ## set rownames to ID column so that ordering of ICs is carried forward
    IC.O <- t(D.O.beta.alt)
    
  }else{ #getSEcoef is NOT called from survMSM (i.e. fit_xMSM functions) or crude survMSM
    
    IC.data <- MSMdata[ ,c(unique(c(ID,clusterID)), t.var, Yname, cum.IPAW, MSMpredict), with = FALSE]
    
    IC.data[, "MSMepsilon" := get(Yname) - get(MSMpredict)]
    
    if (IPW_MSMestimator) {
      sweep.tmp <- IC.data[[cum.IPAW]] * IC.data[["MSMepsilon"]]
    }else{
      sweep.tmp <- IC.data[["MSMepsilon"]]
    }
    
    D.O.beta <- sweep(t(MSMdesign), 2, sweep.tmp, "*") # multiply each column of t(MSMdesign) by sweep.tmp:
    xDT <- data.table(ID = IC.data[[clusterID]], t(D.O.beta)) ## sum of t and theta with no clustering and over t, theta, and cluster observations otherwise
    setkeyv(xDT, cols = "ID")
    xDT <- xDT[, lapply(.SD, sum), by = ID] ## take sum by ID, no longer drop ID column
    cols <- setdiff(names(xDT), "ID")
    xDT[,(cols) := nID/clusterNperID[,sum(clusterN)]*.SD, .SDcols = cols] # apply the clustering weights (n/N where n is iid sample size and N = sum_cluster of cluster size
    D.O.beta.alt <- as.matrix(xDT, rownames = "ID") ## set rownames to ID column so that ordering of ICs is carried forward
    D.O.beta.alt <- t(D.O.beta.alt)
    
    if (gvars$verbose) {
      print("Should all be very small: "); print(abs(apply(D.O.beta.alt, 1, mean)))
    }
    
    if(IPW_MSMestimator & family=="binomial"){
      sweep.tmp <- IC.data[[cum.IPAW]] * IC.data[[MSMpredict]] * (1 - IC.data[[MSMpredict]])
    }else if(IPW_MSMestimator & family=="gaussian"){
      sweep.tmp <- IC.data[[cum.IPAW]]
    }else if(IPW_MSMestimator & family=="poisson"){
      sweep.tmp <- IC.data[[cum.IPAW]] * IC.data[[MSMpredict]]
    }else if(!IPW_MSMestimator & family=="binomial"){
      sweep.tmp <- IC.data[[MSMpredict]] * (1 - IC.data[[MSMpredict]])
    }else{
      stop("IC for coefficients of an MSM of this family is not yet implemented")
    }
    
    C <- sweep(t(MSMdesign), 2, sweep.tmp, "*")
    C <- (1/clusterNperID[,sum(clusterN)]) * (C %*% MSMdesign) # also apply cluster weight nID/clusterNperID[,sum(clusterN)] but numerator cancles out with denominator (not shown) = nID
    Cm1 <- solve(C)
    IC.O <- Cm1 %*% D.O.beta.alt # rownames of influence curves correspond to associated cluster IDs
  }
  
  if (gvars$verbose) {
    print("Should all be very small: "); print(abs(apply(IC.O, 1, mean)))
  }
  
  var.beta <- (IC.O %*% t(IC.O)) / (nID^2)
  se.beta <- sqrt(diag(var.beta))
  return(list(IC.O = IC.O, se.beta = se.beta))
}

##############################################################
## SE of IPW estimator of Sd(t) for a given d and all t's
##############################################################
### IC.O below is the first output of previous function:
getSE.S <- function(nID, S.d.t.predict, h.d.t.predict, design.d.t, IC.O){
  ### SE for S(t) for all t's
  h.by.dl.dt <- matrix(NA, nrow = length(S.d.t.predict), ncol = ncol(design.d.t))
  for(t.val in 1:length(S.d.t.predict)) {
    h.by.dl.dt[t.val,] <- h.d.t.predict[t.val] * design.d.t[t.val, , drop=FALSE]
  }
  
  sum.h.by.dl.dt <- apply(h.by.dl.dt, 2, cumsum)
  S.by.sum.h.by.dl.dt <- sweep(sum.h.by.dl.dt, 1, S.d.t.predict, "*")
  IC.S <- -S.by.sum.h.by.dl.dt %*% IC.O

  if (gvars$verbose) {
    print("Should all be very small: "); print(abs(apply(IC.S, 1, mean)))
  }
  
  var.S <- (IC.S %*% t(IC.S)) / (nID^2)
  se.S <- sqrt(diag(var.S))
  return(list(IC.S = IC.S, se.S = se.S))
}

##############################################################
## SE of IPW estimator of RD=Sd1(t)-Sd2(t) for two given d1 and d2 and all t's
##############################################################
### IC.S.d1 and IC.S.d2 are the first output of previous function above applied to two different d's:
getSE.RD.d1.minus.d2 <- function(nID, IC.S.d1, IC.S.d2){
  IC.RD <- IC.S.d2 - IC.S.d1 # note that we do IC.d2 minus IC.d1 here even though this if for the RD defined as d1 minus d2 ( we flip the order)
  var.RD <- (IC.RD %*% t(IC.RD)) / (nID^2)
  se.RD <- sqrt(diag(var.RD))
  return(se.RD)
}

