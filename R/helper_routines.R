# ---------------------------------------------------------------------------------------------
# (OLD) Helper routine to obtain the follow-up time -> Needs to be re-written
# ---------------------------------------------------------------------------------------------
#' @export
OLD_get.T.tilde <- function(data,IDname,Yname,Cname,tname){
  T.tilde <- unlist(as.list(by(data,data[,IDname],function(x){
    return( min( x[x[,Yname]==1 & x[,tname]>0, tname]-1, x[x[,Cname]==1,tname] ,na.rm=TRUE) )
  })))
  return(T.tilde)
}

# ----------------------------------------------------------------
#### SUMMARY STATISTICS
# ----------------------------------------------------------------
makeFreqTable <- function(rawFreq){
  ntot <- sum(rawFreq)
  rawFreqPercent <- rawFreq/ntot*100
  fineFreq <- cbind("Frequency"=rawFreq,"\\%"=round(rawFreqPercent,2),"Cumulative Frequency"=cumsum(rawFreq),"Cumulative \\%"=round(cumsum(rawFreqPercent),2))
  fineFreq <- as.data.frame(fineFreq)
  eval(parse(text=paste('fineFreq <- cbind(',deparse(substitute(rawFreq)),'=factor(c("',paste(names(rawFreq),collapse='","'),'")),fineFreq)',sep="")))
  rownames(fineFreq) <- 1:nrow(fineFreq)
  fineFreq <- as.matrix(fineFreq)
  return(fineFreq)
}

# ----------------------------------------------------------------
#### SUMMARY STATISTICS
# ----------------------------------------------------------------
makeSumFreqTable <- function(x.freq,cutoffs,varName){
  if(class(cutoffs)=="Date"){
    x.values <- as.Date(names(x.freq))
  }else{
    x.values <- names(x.freq)
  }
  na.yes <- as.logical(sum(is.na( as.numeric(x.values))))
  x.freq.sum <- rep(NA,length(cutoffs)+1+as.numeric(na.yes))
  x.freq.sum[1] <- sum(x.freq[ as.numeric(x.values) < cutoffs[1] ], na.rm=TRUE)
  for(i in 1:(length(cutoffs)-1))
    x.freq.sum[1+i] <- sum(x.freq[ as.numeric(x.values) >= cutoffs[i] & as.numeric(x.values)<cutoffs[i+1] ], na.rm=TRUE)
  x.freq.sum[length(cutoffs)+1] <- sum(x.freq[ as.numeric(x.values)>=cutoffs[length(cutoffs)] ], na.rm=TRUE)
  if(na.yes){
      x.freq.sum[length(cutoffs)+2] <- x.freq[ is.na(as.numeric(x.values)) ]
  }

  catNames <- rep(NA,length(cutoffs)+1+as.numeric(na.yes))
  catNames[1] <- paste("<",cutoffs[1],sep="")
  for(i in 1:(length(cutoffs)-1))
    catNames[i+1] <- paste("[",cutoffs[i],", ",cutoffs[i+1],"[",sep="")
  catNames[length(cutoffs)+1] <- paste(">=",cutoffs[length(cutoffs)],sep="")
  if(na.yes)catNames[length(cutoffs)+2] <- "Missing"
  names(x.freq.sum) <- catNames

  x.freq.sum <- makeFreqTable(x.freq.sum)
  colnames(x.freq.sum)[1] <- varName
  x.freq.sum[length(cutoffs)+1,1] <- paste("$\\geq$ ",cutoffs[length(cutoffs)],sep="")
  return(x.freq.sum)
}


make.table.m0 <- function(S.IPAW, RDscale = "-" , nobs = 0, esti = "IPAW", t.period, se.RDscale.Sdt.K){
  if (missing(se.RDscale.Sdt.K)) {
    se.RDscale.Sdt.K <- matrix(NA, nrow = length(S.IPAW), ncol = length(S.IPAW))
    colnames(se.RDscale.Sdt.K) <- names(S.IPAW)
    rownames(se.RDscale.Sdt.K) <- names(S.IPAW)
  }
  dtheta <- names(S.IPAW)
  RDtable <- matrix(NA, nrow = factorial(length(dtheta)), ncol = length(dtheta)-1)
  # RDtable <- matrix(NA,nrow=2*3,ncol=3)

  ContrastScale <- ifelse(RDscale,"-","/")

  H0val <- ifelse(RDscale,0,1)
  ## Compare mean bootstrap to point estimates - should be similar
  PYK1.IPAW <- allRDtable <- vector("list",length(S.IPAW))
  names(PYK1.IPAW) <- names(allRDtable) <- names(S.IPAW)
  PYK1.IPAW <- lapply(S.IPAW,function(x,y)return(1-x[y+1]),y=t.period)

  rownames(RDtable) <- rep(rev(dtheta)[-length(dtheta)],each=2)
  colnames(RDtable) <- dtheta[-length(dtheta)]

  for(rule1 in rev(dtheta)[-length(dtheta)]) {
    for(rule2 in dtheta[ 1:(which(dtheta==rule1)-1) ]){
      RD <- mapply(ContrastScale,PYK1.IPAW[[rule1]],PYK1.IPAW[[rule2]])
      (pt <- round(RD,4))
      (CI <- round(RD+c(qnorm(0.025),-qnorm(0.025))*se.RDscale.Sdt.K[rule1,rule2],4))
      (ptCI <- paste(pt," [",CI[1],";",CI[2],"]",sep=""))
      (pval <- paste("SE=",round(se.RDscale.Sdt.K[rule1,rule2],4),", p=", round(2*pnorm( abs((RD-H0val)/se.RDscale.Sdt.K[rule1, rule2]), lower=FALSE ),2) ,sep=""))
      RDtable[ rownames(RDtable)%in%rule1,rule2] <- c(ptCI,pval)
    }
  }

  # browser()
  # (seq_along(dtheta)-1) *

  RDtable[is.na(RDtable)] <- ""
  RDtable <- cbind("d1 (row) | d2 (col)"=rep(rev(dtheta)[-length(dtheta)],each=2) , RDtable)
  RDtable[,1] <- paste(gsub("d","$d_{",RDtable[,1]),sep="")
  RDtable[,1] <- paste(gsub("g","}g_{",RDtable[,1]),"}$",sep="")

  # RDtable[ c(0,1,2) * 2 + rep(2, each = 3), 1] <- ""
  # colnames(RDtable)[-1] <- paste(gsub("d","$d_{",colnames(RDtable)[-1]),sep="")
  # colnames(RDtable)[-1] <- paste(gsub("g","}g_{",colnames(RDtable)[-1]),"}$",sep="")
  # colnames(RDtable)[1] <- c("$d_1$ (row) - $d_2$ (col)")

  fileText <- paste(esti,ifelse(RDscale,"RD","RR"),sep="")
  captionText2 <- ifelse(RDscale,"differences","ratios")
  captionText3 <- ifelse(RDscale,"$-$","$/$")
  if(esti=="IPAW")estimates <- "Stabilized inverse weighting"
  if(esti=="IPAWtrunc")estimates <- "Stabilized, truncated inverse weighting"
  if(esti=="crude")estimates <- "Crude"
  model <- "MSM"
  if(esti=="crude")model <- "model"
  est <- "IPAW"
  ##     browser()

  caption <- paste(estimates,
      " estimates of the (cumulative) RD progression risk ",
      captionText2,", $d_1$",
      captionText3,"$d_2$, (with ",gsub("g","",0),
      " grace period(s) allowed) over ",
      (t.period+1)/4," years (",t.period+1," periods). The risk contrasts are derived from a logistic ", model,
      " for the discrete-time hazards fitted based on ", nobs,
      " observations. Variance estimates are derived based on the influence curve of the estimator.",sep="")

  # sink(file.path(res_outf_newdat, paste0(fileText,est,K,"Cm",0,"hICbased.tex")))
  # latex(RDtable,file="",where="!htpb",colheads=colnames(RDtable),
  #   caption = caption,
  #   label= paste(fileText,est,K,"Cm",0,"hbased",sep=""),booktabs=TRUE,rowname=NULL,landscape=TRUE,n.rgroup=rep(2,3))
  # sink()
  return(list(RDtable = RDtable, caption = caption))
}


# ---------------------------------------------------------------------------------------------
#' Helper routine to convert the data into the format data.table required by estmr() function.
#'
#' @param data Input data.table or data.frame.
#' @param ID The name of the unique subject identifier (character, numeric or factor).
#' @param t The name of the time/period variable in \code{data}.
#' @param I The name of the numeric biomarker value which determines the dynamic treatment rule at each time point t.
#' @param imp.I The name of the binary indicator of missingness or imputation for I at time point t and it is used for coding MONITOR(t-1):=1-imp.I(t).
#'  When imp.I(t)=1 it means that the patient was not observed (no office visit) at time-point t and hence no biomarker was measured.
#' @param MONITOR.name The name of the MONITORing variable which will be evaluated by this routine.
#'  This new column MONITOR(t) is the indicator of being monitored (having a doctors visit) at time-point t+1 the indicator of the
#'  imputation (having observed/measured biomarker) at time-point t+1.
#' @param tsinceNis1 The name of the variable that counts number of periods since last monitoring event at t-1.
#'
#' @section Details:
#'
#' Convert the input long format data with the time ordering: (I(t), imp.I(t), C(t), A(t))
#' into the data format required by the estmr() function: (I(t), C(t), A(t), N(t):=1-imp.I(t+1)).
#' N(t) at time-point t is defined as the indicator of being observed (having an office visit) at time point t+1 (next timepoint after t)
#' The very first value of I(t) (at the first time-cycle) is ALWAYS ASSUMED observed/measured (hence I.imp=0 for each first subject-time observation).
#'
#' @section The format of the input and output datasets:
#'
#' The input data.frame data needs to be in long format.
#'
#' The format of the specified columns needs to be as follows.
#'
#' The time ordering of the input data at each t is as follows: (I, imp.I, CENS, TRT)
#'
#' The time ordering of the output data at each t is as follows: (I, CENS, TRT, MONITOR), where MONITOR(t)=1-imp.I(t+1).
#'
#' In output data.table, MONITOR(t-1)=1 indicates that the biomarker I(t) at t is observed and vice versa.
#'
#' In addition the output data.table will contain a column "tsinceNis1", where:
#' \itemize{
#'   \item tsinceNis1(t) = 0 means that the person was monitored at time-point t-1.
#'   \item tsinceNis1(t) > 0 is the count of the number of cycles since last monitoring event.
#' }
#' @return A data.table in long format with ordering (I, CENS, TRT, MONITOR)
#' @export
convertdata <- function(data, ID, t, imp.I, MONITOR.name = "N", tsinceNis1 = "tsinceNis1"){
  ID.expression <- as.name(ID)
  indx <- as.name("indx")
  if (is.data.table(data)) {
    DT <- data.table(data, key=c(ID, t))
    # DT <- data.table(data[,c(ID, t, imp.I), with = FALSE], key=c(ID, t))
  } else if (is.data.frame(data)) {
    DT <- data.table(data, key=c(ID, t))
    # DT <- data.table(data[,c(ID, t, imp.I)], key=c(ID, t))
  } else {
    stop("input data must be either a data.table or a data.frame")
  }
  # "Leading" (shifting up) and inverting indicator of observing I, renaming it to MONITOR value;
  # N(t-1)=1 indicates that I(t) is observed. Note that the very first I(t) is assumed to be always observed.
  DT[, (MONITOR.name) := shift(.SD, n=1L, fill=NA, type="lead"), by = eval(ID.expression), .SDcols = (imp.I)]
  DT[, (MONITOR.name) := 1L - get(MONITOR.name)]
  # Create "indx" vector that goes up by 1 every time MONITOR.name(t-1) shifts from 1 to 0 or from 0 to 1
  DT[, ("indx") := cumsum(c(FALSE, get(MONITOR.name)!=0L))[-.N], by = eval(ID.expression)]
  DT[, (tsinceNis1) := seq(.N)-1, by = list(eval(ID.expression), eval(indx))]
  # DT[, (tsinceNis1) := seq(.N)-1, by = .(eval(ID.expression), indx)]
  DT[is.na(DT[["indx"]]), (tsinceNis1) := NA]
  DT[, ("indx") := NULL]
  # DT[, (imp.I) := NULL]
  return(DT)
}

# ---------------------------------------------------------------------------------------------
#' Define the indicators of following/not-following specific dynamic treatment rules indexed by theta.
#'
#' @param data Input data.frame or data.table in long format, see below for the description of the assumed format.
#' @param theta The vector of continuous cutoff values that index each dynamic treatment rule
#' @param ID The name of the unique subject identifier
#' @param t The name of the variable indicating time-period
#' @param I Continuous biomarker variable used for determining the treatmet decision rule
#' @param CENS Binary indicator of being censored at t;
#' @param TRT Binary indicator of the treatment (exposure) at t;
#' @param MONITOR The indicator of having a visit and having measured or observed biomarker I(t+1) (the biomarker value at THE NEXT TIME CYCLE).
#'  In other words the value of MONITOR(t-1) (at t-1) being 1 indicates that I(t) at time point t was observed/measured.
#'  The very first value of I(t) (at the first time-cycle) is ALWAYS ASSUMED observed/measured.
#' @param rule.names Vector of column names for indicators of following/not following each rule (must be the same dimension as theta).
#'  When not supplied the following convention is adopted for naming these columns: paste0("d",theta).
#' @param return.allcolumns Set to \code{TRUE} to return the original data columns along with new columns that define each rule
#' (can be useful when employing piping/sequencing operators).
#'
#' @section Details:
#'
#' * This function takes an input data.frame or data.table data
#'   and produces an output data.table with indicators/probabilities of following not following a specific treatment rule.
#'   The resulting data.table is used internally by the stremr() function to determine which observation is following each specific rule at each time point t.
#'
#' * Evaluates which observations were following the dynamic-decision treatment rule defined by the measured biomarker I
#'   and pre-defined cutoffs of the input vector theta.
#'
#' * Produces a separate rule indicator column for each value in the input vector theta based on the following dynamic rule at t:
#'\itemize{
#' \item (1) Follow rule at t if uncensored (C(t)=0) and remaining on treatment (A(t-1)=A(t)=1)
#' \item (2) Follow rule at t if uncensored (C(t)=0), haven't changed treatment and wasn't monitored (MONITOR(t-1)=0)
#' \item (3) Follow rule at t if uncensored (C(t)=0), was monitored (MONITOR(t-1)=0) and either of the two:
#'     (A) (I(t) >= d.theta) and switched to treatment at t; or (B) (I(t) < d.theta) and haven't changed treatment
#'}
#'
#' * The format (time-ordering) of data is the same as required by the stremr() function: (I(t), CENS(t), TRT(t), MONITOR(t)).
#'   MONITOR(t) at time-point t is defined as the indicator of being observed (having an office visit) at time point t+1 (next timepoint after t)
#'   It is assumed that imp.I(t) is always 0 for the very first time-point.
#'
#' @examples
#'
#' \dontrun{
#' theta <- seq(7,8.5,by=0.5)
#' FOLLOW.D.DT <- follow.rule.d.DT(data = data, theta = theta,
#'                  ID = "StudyID", t = "X_intnum", I = "X_a1c",
#'                  TRT = "X_exposure", CENS = "X_censor", MONITOR = "N.t",
#'                  rule.names = paste0("new.d",theta))
#' }
#' @return A data.table with a separate column for each value in \code{theta}. Each column consists of indicators of following/not-following
#'  each rule indexed by a value form \code{theta}. In addition, the returned data.table contains \code{ID} and \code{t} columns for easy merging
#'  with the original data.
#' @export
#'
follow.rule.d.DT <- function(data, theta, ID, t, I, CENS, TRT, MONITOR, rule.names = NULL, return.allcolumns = FALSE){
  # require('data.table')
  ID.expression <- as.name(ID)
  indx <- as.name("indx")
  chgTRT <- as.name("chgTRT")
  lastN.t <- as.name("lastN.t")

  if (return.allcolumns) {
    DT <- data.table(data, key=c(ID,t))
  } else if (is.data.table(data)) {
    DT <- data.table(data[,c(ID, t, TRT, CENS, I, MONITOR), with = FALSE], key=c(ID,t))
  } else if (is.data.frame(data)) {
    DT <- data.table(data[,c(ID, t, TRT, CENS, I, MONITOR)], key=c(ID,t))
  } else {
    stop("input data must be either a data.table or a data.frame")
  }
  # Create "indx" vector that goes up by 1 every time MONITOR(t-1) shifts from 1 to 0 or from 0 to 1
  DT[, "indx" := cumsum(c(FALSE, get(MONITOR)!=0L))[-.N], by = eval(ID.expression)]
  # Intermediate variable lastN.t to count number of cycles since last visit at t-1. Reset lastN.t(t)=0 when MONITOR(t-1)=1.
  DT[, "lastN.t" := seq(.N)-1, by = list(eval(ID.expression), eval(indx))]
  DT[is.na(DT[["indx"]]), "lastN.t" := NA]
  DT[, "indx" := NULL]

  # Define chgTRT=TRT(t)-TRT(t-1), switching to treatment (+1), not changing treatment (0), going off treatment (-1). Assume people were off treatment prior to t=0.
  DT[, "chgTRT" := diff(c(0L, .SD[[1]])), by = eval(ID.expression), .SDcols=(TRT)]
  # The observation is not censored at time-point t
  # DT[, notC := (get(CENS)==0L), by = eval(ID.expression)]

  # (1) Follow rule at t if uncensored and remaining on treatment (TRT(t-1)=TRT(t)=1):
  DT[, "d.follow_r1" := (get(CENS)==0L) & (get(TRT)==1L) & (eval(chgTRT)==0L), by = eval(ID.expression), with = FALSE]
  # rule1: (C[t] == 0L) & (A[t-1] == 1L) & (A[t] == 1)

  # (2) Follow rule at t if uncensored, haven't changed treatment and wasn't monitored (MONITOR(t-1)=0)
  DT[, "d.follow_r2" := (get(CENS)==0L) & (eval(lastN.t) > 0L) & (eval(chgTRT)==0L), by = eval(ID.expression), with = FALSE]
  # rule2: (C[t] == 0L) & (N[t-1] == 0) & (A[t-1] == A[t])

  # (3) Follow rule at t if uncensored, was monitored (MONITOR(t-1)=0) and either:
  # (A) (I(t) >= d.theta) and switched to treatment at t; or (B) (I(t) < d.theta) and haven't changed treatment
  for (dtheta in theta) {
    DT[, "d.follow_r3" := (get(CENS)==0L) & (eval(lastN.t) == 0L) & (((get(I) >= eval(dtheta)) & (eval(chgTRT)==1L)) | ((get(I) < eval(dtheta)) & (eval(chgTRT)==0L))), by = eval(ID.expression), with = FALSE]
    # rule3: (C[t] == 0L) & (N[t-1] == 1) & ((I[t] >= d.theta & A[t] == 1L & A[t-1] == 0L) | (I[t] < d.theta & A[t] == 0L & A[t-1] == 0L))
    # ONE INDICATOR IF FOLLOWING ANY OF THE 3 ABOVE RULES AT each t:
    DT[, "d.follow_allr" := eval(parse(text="d.follow_r1 | d.follow_r2 | d.follow_r3")), by = eval(ID.expression)]
    # INDICATOR OF CONTINUOUS (UNINTERRUPTED) RULE FOLLOWING from t=0 to EOF:
    DT[, paste0("d",dtheta) := eval(parse(text="as.logical(cumprod(d.follow_allr))")), by = eval(ID.expression)]
  }

  DT[, "chgTRT" := NULL]; DT[, "d.follow_r1" := NULL]; DT[, "d.follow_r2" := NULL]; DT[, "d.follow_r3" := NULL]; DT[, "d.follow_allr" := NULL]
  if (!is.null(rule.names)) {
    stopifnot(length(rule.names)==length(theta))
    setnames(DT, old = paste0("d",theta), new = rule.names)
  } else {
    rule.names <- paste0("d",theta)
  }

  if (!return.allcolumns) {
    DT <- DT[, c(ID, t, rule.names), with=FALSE]
  }
  return(DT)
}
