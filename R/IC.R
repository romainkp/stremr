###################################
## SE of IPW estimator of MSM coef
###################################
getSEcoef <- function(ID, nID, t.var, Yname, MSMdata, MSMpredict, MSMdesign, IPW_MSMestimator = TRUE){
  # browser()
  # head(MSMdata)

  cumm.IPAW <- "cumm.IPAW"

  IC.data <- MSMdata[ ,c(ID, t.var, Yname, cumm.IPAW, MSMpredict), with = FALSE] # contains data for all obs compatible with at least one rule
  # IC.data <- MSMdata[ ,c(ID,t.var,Yname,cumm.IPAW), with = FALSE] # contains data for all obs compatible with at least one rule

  # MSMepsilon <- IC.data[,Yname, with = FALSE]-IC.data[,MSMpredict, with = FALSE]
  # IC.data <- cbind(IC.data,MSMepsilon)

  IC.data[, "MSMepsilon" := get(Yname) - get(MSMpredict)]

  if(IPW_MSMestimator){
    sweep.tmp <- IC.data[[cumm.IPAW]] * IC.data[["MSMepsilon"]]
  }else{
    sweep.tmp <- IC.data[, "MSMepsilon"]
  }
  # *********
  #OS: multiply each column of t(MSMdesign) by sweep.tmp:
  D.O.beta <- sweep(t(MSMdesign), 2, sweep.tmp, "*")

  xDT <- data.table(ID = IC.data[[ID]], t(D.O.beta))
  setkeyv(xDT, cols = "ID")
  D.O.beta.alt <- as.matrix(xDT[, lapply(.SD, sum), by = ID][, ID := NULL])
  D.O.beta.alt <- t(D.O.beta.alt)

  # old approach:
  # t.eval <- system.time(
  #   D.O.beta <- apply(D.O.beta, 1, function(x, y) {
  #     return(tapply(x, y, sum))
  #     # browser()
  #     # res1 <- tapply(x, y, sum)
  #     # res2 <- IC.data[,sum(x), by = ID][["V1"]]
  #     # res1[1]
  #     # res2[2]
  #   }, y = as.factor(IC.data[[ID]]))
  # )
  # D.O.beta <- t(D.O.beta)

  # head(D.O.beta)
  # dim(D.O.beta)
  # print("t.eval tapply: "); print(t.eval)
  # for 10K:
  #  user  system elapsed
  # 0.518   0.050   0.583

  print("Should all be very small: "); print(abs(apply(D.O.beta.alt, 1, mean)))
  # print("Should all be very small: "); print(abs(apply(D.O.beta, 1, mean)))

  if(IPW_MSMestimator){
    sweep.tmp <- IC.data[[cumm.IPAW]] * IC.data[[MSMpredict]] * (1-IC.data[[MSMpredict]])
  } else {
    sweep.tmp <- IC.data[[MSMpredict]] * (1-IC.data[[MSMpredict]])
  }

  C <- sweep(t(MSMdesign), 2, sweep.tmp, "*")
  C <- (C %*% MSMdesign) / nID
  Cm1 <- solve(C)

  IC.O <- Cm1 %*% D.O.beta.alt
  # IC.O <- Cm1 %*% D.O.beta

  print("Should all be very small: "); print(abs(apply(IC.O, 1, mean)))

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

  print("Should all be very small: "); print(abs(apply(IC.S, 1, mean)))
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

