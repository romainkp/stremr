###################################
## SE of IPW estimator of MSM coef
###################################
getSEcoef <- function(ID,t.var,MSMdata,MSMpredict,MSMdesign,Yname,MSMestimator,nID){
  IC.data <- MSMdata[ ,c(ID,t.var,Yname,MSMestimator)] # contains data for all obs compatible with at least one rule
  IC.data <- cbind(IC.data,MSMpredict)
  MSMepsilon <- IC.data[,Yname]-IC.data[,"MSMpredict"]
  IC.data <- cbind(IC.data,MSMepsilon)
  if(MSMestimator%in%c("IPAWtrunc","IPAW")){
    sweep.tmp <- IC.data[,MSMestimator]*IC.data[,"MSMepsilon"]
  }else{
    sweep.tmp <- IC.data[,"MSMepsilon"]
  }
  # *********
  #OS: multiply each column of t(MSMdesign) by sweep.tmp:
  D.O.beta <- sweep(t(MSMdesign),2,sweep.tmp,"*")
  D.O.beta <- apply(D.O.beta,1,function(x,y){
    return(tapply(x,y,sum))
  },y=as.factor(IC.data[,ID]))
  D.O.beta <- t(D.O.beta)
  cat("Should all be very small: ",abs(apply(D.O.beta,1,mean)))

  if(MSMestimator%in%c("IPAWtrunc","IPAW")){
    sweep.tmp <- IC.data[,MSMestimator]*IC.data[,"MSMpredict"]*(1-IC.data[,"MSMpredict"])
  }else{
    sweep.tmp <- IC.data[,"MSMpredict"]*(1-IC.data[,"MSMpredict"])
  }
  C <- sweep(t(MSMdesign),2, sweep.tmp ,"*")
  C <- (C%*%MSMdesign)/nID
  Cm1 <- solve(C)

  IC.O <- Cm1%*%D.O.beta
  cat("Should all be very small: ",abs(apply(IC.O,1,mean)))
  var.beta <- (IC.O%*%t(IC.O))/(nID^2)
  se.beta <- sqrt(diag(var.beta))
  return(list(IC.O,se.beta))
}

##############################################################
## SE of IPW estimator of Sd(t) for a given d and all t's
##############################################################
### IC.O below is the first output of previous function:
getSE.S <- function(nID,S.d.t.predict,h.d.t.predict,design.d.t,IC.O){
  ### SE for S(t) for all t's
  h.by.dl.dt <- matrix(NA,nrow=length(S.d.t.predict),ncol=ncol(design.d.t))
  for(t.val in 1:length(S.d.t.predict))
    h.by.dl.dt[t.val,] <- h.d.t.predict[t.val]*design.d.t[t.val,,drop=FALSE]
  sum.h.by.dl.dt <- apply(h.by.dl.dt,2,cumsum)
  S.by.sum.h.by.dl.dt <- sweep(sum.h.by.dl.dt,1,S.d.t.predict,"*")
  IC.S <- -S.by.sum.h.by.dl.dt%*%IC.O
  cat("Should all be very small: ",abs(apply(IC.S,1,mean)))
  var.S <- (IC.S%*%t(IC.S))/(nID^2)
  se.S <- sqrt(diag(var.S))
  return(list(IC.S,se.S))
}

##############################################################
## SE of IPW estimator of RD=Sd1(t)-Sd2(t) for two given d1 and d2 and all t's
##############################################################
### IC.S.d1 and IC.S.d2 are the first output of previous function above applied to two different d's:
getSE.RD.d1.minus.d2 <- function(nID,IC.S.d1,IC.S.d2){
  IC.RD <- IC.S.d2 - IC.S.d1 # note that we do IC.d2 minus IC.d1 here even though this if for the RD defined as d1 minus d2 ( we flip the order)
  var.RD <- (IC.RD%*%t(IC.RD))/(nID^2)
  se.RD <- sqrt(diag(var.RD))
  return(se.RD)
}

