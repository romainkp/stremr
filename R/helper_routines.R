# ---------------------------------------------------------------------------------------------
# (OLD) Helper routine to obtain the follow-up time -> Needs to be re-written
# ---------------------------------------------------------------------------------------------
OLD_get.T.tilde <- function(O.data,IDname,Yname,Cname,tname){
  T.tilde <- unlist(as.list(by(O.data,O.data[,IDname],function(x){
    return( min( x[x[,Yname]==1 & x[,tname]>0, tname]-1, x[x[,Cname]==1,tname] ,na.rm=TRUE) )
  })))
  return(T.tilde)
}

# ---------------------------------------------------------------------------------------------
# Helper routine to convert the data into the format data.table required by estmr() function.
# ---------------------------------------------------------------------------------------------
# (1) Convert the input long format data with the time ordering: (I(t), imp.I(t), C(t), A(t))
#     into the data format required by the estmr() function: (I(t), C(t), A(t), N(t):=1-imp.I(t+1)).
#     N(t) at time-point t is defined as the indicator of being observed (having an office visit) at time point t+1 (next timepoint after t)
#     The very first value of I(t) (at the first time-cycle) is ALWAYS ASSUMED observed/measured (hence I.imp=0 for each first subject-time observation).
# (2) Define a new covariate lastN.t(t) as part of L(t), which counts the number of time-points (cycles) since last visit at t-1:
#     * lastN.t(t) = 0 means that the person was monitored at time-point t-1.
#     * lastN.t(t) > 0 is the count of the number of cycles since last monitoring event.
# ---------------------------------------------------------------------------------------------
# The format of the input and output datasets.
# ---------------------------------------------------------------------------------------------
# The input data.frame Odata needs to be in long format.
# The format of the specified columns needs to be as follows.
# The time ordering of the input data at each t is as follows: (I, imp.I, C, A)
# The time ordering of the output data at each t is as follows: (I, C, A, N), where N(t)=1-imp.I(t+1).
# In output data.table, N(t-1)=1 indicates that the biomarker I(t) at t is observed and vice versa.
# In addition the output data.table will contain a column "lastN.t", where:
#   * lastN.t(t) = 0 means that the person was monitored at time-point t-1.
#   * lastN.t(t) > 0 is the count of the number of cycles since last monitoring event.
# ---------------------------------------------------------------------------------------------
# Arguments:
# ---------------------------------------------------------------------------------------------
# Odata - input data.table or data.frame
# ID - the name of the unique subject identifier (character, numeric or factor)
# I - the name of the numeric biomarker value which determines the dynamic treatment rule at each time point t.
# imp.I - the name of the binary indicator of missingness or imputation for I at time point t and it is used for for coding N(t-1):=1-imp.I(t).
#     When imp.I(t)=1 it means that the patient was not observed (no office visit) at time-point t and hence no biomarker was measured.
# N.name - the name of the variable which will be evaluated by this routine.
#     This new column MONITOR(t) is the indicator of being monitored (having a doctors visit) at time-point t+1 the indicator of the
#     imputation (having observed/measured biomarker) at time-point t+1
convertOdata <- function(Odata, ID, t, imp.I, N.name){
  require('data.table')
  if (is.data.table(Odata)) {
    DT <- data.table(Odata[,c(ID, t, imp.I), with = FALSE], key=c(ID, t))
  } else if (is.data.frame(Odata)) {
    DT <- data.table(Odata[,c(ID, t, imp.I)], key=c(ID, t))
  } else {
    stop("input Odata must be either a data.table or a data.frame")
  }
  # "Leading" (shifting up) and inverting indicator of observing I, renaming it to MONITOR value;
  # N(t-1)=1 indicates that I(t) is observed. Note that the very first I(t) is assumed to be always observed.
  DT[, (N.name) := shift(.SD, n=1L, fill=NA, type="lead"), by=get(ID), .SDcols=(imp.I)]
  DT[, (N.name) := 1L - get(N.name)]
  # Create "indx" vector that goes up by 1 every time N.name(t-1) shifts from 1 to 0 or from 0 to 1
  DT[, indx:=cumsum(c(FALSE, get(N.name)!=0L))[-.N], by = get(ID)]
  DT[, lastN.t:=seq(.N)-1, by = .(Study_ID, indx)]
  DT[is.na(DT[["indx"]]), lastN.t:=NA]
  DT[, indx:=NULL]
  return(DT)
}

# ---------------------------------------------------------------------------------------------
# Define the indicators of following/not-following specific dynamic treatment rules indexed by theta.
# ---------------------------------------------------------------------------------------------
# * This function takes an input data.frame or data.table Odata
#   and produces an output data.table with indicators/probabilities of following not following a specific treatment rule.
#   The resulting data.table is used internally by the estimtr() function to determine which observation is following each specific rule at each time point t.
# * Evaluates which observations were following the dynamic-decision treatment rule defined by the measured biomarker I
#   and pre-defined cutoffs of the input vector theta.
# * Produces a separate rule indicator column for each value in the input vector theta based on the following dynamic rule at t:
# (1) Follow rule at t if uncensored (C(t)=0) and remaining on treatment (A(t-1)=A(t)=1)
# (2) Follow rule at t if uncensored (C(t)=0), haven't changed treatment and wasn't monitored (MONITOR(t-1)=0)
# (3) Follow rule at t if uncensored (C(t)=0), was monitored (MONITOR(t-1)=0) and either of the two:
#     (A) (I(t) >= d.theta) and switched to treatment at t; or (B) (I(t) < d.theta) and haven't changed treatment
# * The format (time-ordering) of Odata is the same as required by the estimtr() function: (I(t), CENS(t), TRT(t), MONITOR(t)).
#   MONITOR(t) at time-point t is defined as the indicator of being observed (having an office visit) at time point t+1 (next timepoint after t)
#   It is assumed that imp.I(t) is always 0 for the very first time-point.
# ---------------------------------------------------------------------------------------------
# Arguments:
# ---------------------------------------------------------------------------------------------
# theta - the vector of continuous cutoff values that index each dynamic treatment rule
# Odata - input data.frame in long format, see below for the description of the assumed format.
# ID - the name of the unique subject identifier
# t - the name of the variable indicating time-period
# I - continuous biomarker variable used for determining the treatmet decision rule
# CENS - binary indicator of being censored at t;
# TRT - binary indicator of the treatment (exposure) at t;
# MONITOR - The indicator of having a visit and having measured or observed biomarker I(t+1) (the biomarker value at THE NEXT TIME CYCLE).
#     In other words the value of MONITOR(t-1) (at t-1) being 1 indicates that I(t) at time point t was observed/measured.
#     The very first value of I(t) (at the first time-cycle) is ALWAYS ASSUMED observed/measured.
# rule.names - vector of column names for indicators of following/not following each rule (must be the same dimension as theta).
#   When not supplied the following convention is adopted for naming these columns: paste0("d",theta).
# ---------------------------------------------------------------------------------------------
# Usage:
# ---------------------------------------------------------------------------------------------
# theta <- seq(7,8.5,by=0.5)
# t.run <- system.time(
#   FOLLOW.D.DT <- follow.rule.d.DT(theta = theta, Odata = Odata,
#                  ID = "Study_ID", t = "X_intnum",  I = "X_a1c",
#                  TRT = "X_exposure", CENS = "X_censor", MONITOR = "N.t",
#                  rule.names = paste0("new.d",theta))
#   )
follow.rule.d.DT <- function(theta, Odata, ID, t, I, CENS, TRT, MONITOR, rule.names = NULL){
  require('data.table')
  if (is.data.table(Odata)) {
    DT <- data.table(Odata[,c(ID, t, TRT, CENS, I, MONITOR), with = FALSE], key=c(ID,t))
  } else if (is.data.frame(Odata)) {
    DT <- data.table(Odata[,c(ID, t, TRT, CENS, I, MONITOR)], key=c(ID,t))
  } else {
    stop("input Odata must be either a data.table or a data.frame")
  }
  # Create "indx" vector that goes up by 1 every time MONITOR(t-1) shifts from 1 to 0 or from 0 to 1
  DT[, indx:=cumsum(c(FALSE, get(MONITOR)!=0L))[-.N], by = get(ID)]
  # Intermediate variable lastN.t to count number of cycles since last visit at t-1. Reset lastN.t(t)=0 when MONITOR(t-1)=1.
  DT[, lastN.t:=seq(.N)-1, by = .(Study_ID, indx)]
  DT[is.na(DT[["indx"]]), lastN.t:=NA]
  DT[, indx:=NULL]

  # Define chgTRT=TRT(t)-TRT(t-1), switching to treatment (+1), not changing treatment (0), going off treatment (-1). Assume people were off treatment prior to t=0.
  DT[, "chgTRT":=diff(c(0L, .SD[[1]])), by = get(ID), .SDcols=(TRT)]
  # The observation is not censored at time-point t
  DT[, notC:=(get(CENS)==0L), by = get(ID)]

  # (1) Follow rule at t if uncensored and remaining on treatment (TRT(t-1)=TRT(t)=1):
  DT[, "d.follow_r1":= notC & (get(TRT)==1L) & (chgTRT==0L), by = get(ID), with = FALSE]
  # (2) Follow rule at t if uncensored, haven't changed treatment and wasn't monitored (MONITOR(t-1)=0)
  DT[, "d.follow_r2":= notC & (lastN.t > 0L) & (chgTRT==0L), by = get(ID), with = FALSE]
  # (3) Follow rule at t if uncensored, was monitored (MONITOR(t-1)=0) and either:
  # (A) (I(t) >= d.theta) and switched to treatment at t; or (B) (I(t) < d.theta) and haven't changed treatment
  for (dtheta in theta) {
    DT[, "d.follow_r3":= notC & (lastN.t == 0L) & (((get(I) >= eval(dtheta)) & (chgTRT==1L)) | ((get(I) < eval(dtheta)) & (chgTRT==0L))), by = get(ID), with = FALSE]
    # ONE INDICATOR IF FOLLOWING ANY OF THE 3 ABOVE RULES AT each t:
    DT[, "d.follow_allr":= d.follow_r1 | d.follow_r2 | d.follow_r3, by = get(ID)]
    # INDICATOR OF CONTINUOUS (UNINTERRUPTED) RULE FOLLOWING from t=0 to EOF:
    DT[, paste0("d",dtheta):=as.logical(cumprod(d.follow_allr)), by = get(ID)]
  }

  DT[, d.follow_r3:=NULL]; DT[, d.follow_allr:=NULL]
  if (!is.null(rule.names)) {
    stopifnot(length(rule.names)==length(theta))
    setnames(DT, old = paste0("d",theta), new = rule.names)
  } else {
    rule.names <- paste0("d",theta)
  }
  return(DT[, c(ID, t, rule.names), with=FALSE])
}
