# # --------------------------------
# # Define event indicator and right-censoring indicator, define total follow-up length
# # --------------------------------
# define_indicators <- function(O.data, ID = "ID", t = "t", TRT = "TI", CENS = "C", MONITOR = "N", I = "highA1c") {
#   # Add indicator Delta=I(T<=C) (Y is 1 at some point for the subject):
#   O.dataDT <- data.table(O.data, key = c(ID, t))
#   EVENT_IND <- "Delta";
#   O.dataDT[,(EVENT_IND) := as.integer(any(get(outcome) %in% 1)), by=eval(ID)]
#   # Add indicator Delta=I(T>C) (subject was right-censored at some point):
#   CENS_IND <- "AnyCensored"; noCENScat <- 0L; CENS <- c("C")
#   O.dataDT[, (CENS_IND) := FALSE, by=eval(ID)]
#   for (Cvar in CENS) {
#     O.dataDT[, (CENS_IND) := get(CENS_IND) | any(!get(Cvar) %in% c(eval(noCENScat),NA)), by = eval(ID)]
#   }
#   O.dataDT[, (CENS_IND) := as.integer(get(CENS_IND))]
#   ## Add variable that indicates if TI initiation occured previously, relevant for model of TI continuation
#   TIcovarname <- "barTIm1eq0"
#   O.dataDT[, (TIcovarname) := as.integer(c(0, cumsum(get(TRT))[-.N]) %in% 0), by = eval(ID)]
#   # O.dataDT[1:100,]
#   return(O.dataDT)
# }
# # ---------------------------------------------------------------------------------------------------------
# # Get indicators of rule followers for two rules (dlow, dhigh) based on binary I:
# # ---------------------------------------------------------------------------------------------------------
# add.gtheta.to.O <- function(O.data, ID = "ID", t = "t", TRT = "TI", CENS = "C", MONITOR = "N", I = "highA1c", rule_name1 = "dlow", rule_name2 = "dhigh"){
#   ## Count how many patients follow each rule; follow rule dlow - start TI at low A1c threshold, i.e. start at study entry.
#   O.dataDT <- data.table(O.data, key = c(ID, t))
#   # Add indicator for rule (dlow) -> start TI when highA1c>=0 and stay on treatment ALL the time (i.e., start at study entry)
#   # equivalent to counterfactual regiment \barA=1
#   rule_name1 <- "dlow"
#   O.dataDT[,(rule_name1) := as.integer(all(get(TRT) %in% c(1,NA))), by=eval(ID)]
#   # set last row to 0, since its always either Y=1 or C=1
#   O.dataDT[O.dataDT[,.I[.N], by = eval(ID)][["V1"]],(rule_name1) := 0L]
#   # set All NAs for treatment (TRT) to be NA for the rule as well
#   O.dataDT[O.dataDT[,.I[which(is.na(get(TRT)))], by = eval(ID)][["V1"]],(rule_name1) := NA]

#   # counterfactual TRT for rule dlow (equivalent to always treated):
#   O.dataDT[,"TI.gstar." %+% rule_name1 := 1L]

#   # Add indicator for rule (dhigh) -> start TRT only when I=1 (highA1c)
#   rule_name2 <- "dhigh"
#   O.dataDT_dhigh <- stremr::defineTRTrules(O.dataDT, theta = 1, ID = ID,
#                                               t = t, I = I, CENS = CENS, TRT = TRT,
#                                               MONITOR = MONITOR,
#                                               tsinceNis1 = "lastNat1",
#                                               rule.names = rule_name2)
#   O.dataDT_dhigh[, (rule_name2) := as.integer(get(rule_name2))]
#   O.dataDT <- merge(O.dataDT, O.dataDT_dhigh, by=c(ID, t))

#   O.dataDT_TIdhigh <- stremr::defineIntervedTRT(O.dataDT, theta = 1, ID = ID,
#                                             t = t, I = I, CENS = CENS, TRT = TRT,
#                                             MONITOR = MONITOR,
#                                             tsinceNis1 = "lastNat1",
#                                             new.TRT.names = "TI.gstar." %+% rule_name2)
#   O.dataDT <- merge(O.dataDT, O.dataDT_TIdhigh, by=c(ID, t))
#   O.dataDT[1:100, ]
#   return(O.dataDT)
# }
# --------------------------------
# (II) Define event indicator and right-censoring indicator, define total follow-up length
# --------------------------------
# O.dataDT <- define_indicators(O.data)
# --------------------------------
# (III) Define rules for dlow & dhigh
# --------------------------------
# O.dataDTrules <- add.gtheta.to.O(O.dataDT)


# -------------------------------------------------------------------------------------------
# Shift the outcome up by 1 and drop all observations that follow afterwards (all NA)
# -------------------------------------------------------------------------------------------
OUTCOME <- "Y"
shifted.OUTCOME <- OUTCOME%+%".tplus1"
O.dataDTrules_Nstar[, (shifted.OUTCOME) := shift(get(OUTCOME), n = 1L, type = "lead"), by = ID]
O.dataDTrules_Nstar <- O.dataDTrules_Nstar[!get(OUTCOME)%in%1,]
# O.dataDTrules_Nstar[1:50]
