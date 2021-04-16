#-------------------------------------------------------------------
# EXAMPLE BASED ON SIMULATED DATA
#-------------------------------------------------------------------
require("data.table")
require("magrittr")
data(OdataCatCENS)
OdataDT <- as.data.table(OdataCatCENS, key=c("ID", "t"))

#-------------------------------------------------------------------
# Define the counterfactual dynamic treatment assignment
#-------------------------------------------------------------------
# Define two dynamic rules: dlow & dhigh
OdataDT <- defineIntervedTRT(OdataDT, theta = c(0,1), ID = "ID", t = "t", I = "highA1c",
                          CENS = "C", TRT = "TI", MONITOR = "N", tsinceNis1 = "lastNat1",
                          new.TRT.names = c("dlow", "dhigh"), return.allcolumns = TRUE)
