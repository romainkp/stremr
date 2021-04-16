context("readme.md example")

options(stremr.verbose = FALSE)
options(gridisl.verbose = FALSE)
options(sl3.verbose = FALSE)
options(condensier.verbose = FALSE)
# options(stremr.verbose = TRUE)
# options(gridisl.verbose = TRUE)
# options(sl3.verbose = TRUE)
# options(condensier.verbose = TRUE)

require("stremr")
require("data.table")
require("magrittr")
data(OdataNoCENS)
OdataDT <- as.data.table(OdataNoCENS, key=c("ID", "t"))

OdataDT[, ("N.tminus1") := shift(get("N"), n = 1L, type = "lag", fill = 1L), by = ID]
OdataDT[, ("TI.tminus1") := shift(get("TI"), n = 1L, type = "lag", fill = 1L), by = ID]

OdataDT[, ("TI.set1") := 1L]
OdataDT[, ("TI.set0") := 0L]

OData <- importData(OdataDT, ID = "ID", t = "t", covars = c("highA1c", "lastNat1", "N.tminus1"), CENS = "C", TRT = "TI", OUTCOME = "Y.tplus1")

##extract data from the data object with helper fun and set the vars
get_data(OData)[, ("TI.set0") := 1L]
get_data(OData)[, ("TI.set0") := 0L]

test_that("get_data can be used to set new column(s) in OData$dat.sVar", {
  expect_equal(unique(OData$dat.sVar[["TI.set0"]]), 0L)
})

gform_CENS <- "C ~ highA1c + lastNat1"
gform_TRT <- "TI ~ CVD + highA1c + N.tminus1"

stratify_CENS <- list(C=c("t < 16", "t == 16"))

test_that("readme examples run as expected without errors", {
  print(OData)
  OData <- fitPropensity(OData,
                         gform_CENS = gform_CENS,
                         gform_TRT = gform_TRT,
                         stratify_CENS = stratify_CENS)

  AKME.St.1 <- getIPWeights(OData, intervened_TRT = "TI.set1") %>%
               survNPMSM(OData) %$%
               estimates

  AKME.St.1[]

  IPW.St.1 <- getIPWeights(OData, intervened_TRT = "TI.set1") %>%
              directIPW(OData) %$%
              estimates
  IPW.St.1[]

  wts.DT.1 <- getIPWeights(OData = OData, intervened_TRT = "TI.set1", rule_name = "TI1")
  wts.DT.0 <- getIPWeights(OData = OData, intervened_TRT = "TI.set0", rule_name = "TI0")
  survMSM_res <- survMSM(list(wts.DT.1, wts.DT.0), OData, tbreaks = c(1:8,12,16)-1,)

  survMSM_res[["TI0"]][["estimates"]]
  survMSM_res[["TI1"]][["estimates"]]

  t.surv <- c(0:15)
  Qforms <- rep.int("Qkplus1 ~ CVD + highA1c + N + lastNat1 + TI + TI.tminus1", (max(t.surv)+1))

  params = defModel(estimator = "speedglm__glm")

  gcomp_est <- fit_GCOMP(OData, tvals = t.surv, intervened_TRT = "TI.set1", Qforms = Qforms, models = params, stratifyQ_by_rule = FALSE)
  gcomp_est$estimates[]

  tmle_est <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "TI.set1", Qforms = Qforms, models = params, stratifyQ_by_rule = TRUE)
  tmle_est$estimates[]

  # require("doParallel")
  # registerDoParallel(cores = parallel::detectCores())
  # tmle_est <- fit_TMLE(OData, tvals = t.surv, intervened_TRT = "TI.set1", Qforms = Qforms, models = params, stratifyQ_by_rule = TRUE, parallel = TRUE)

  # if ("package:h2o" %in% search()) { detach("package:h2o", unload=TRUE) }
  # if ("h2o" %in% rownames(installed.packages())) { remove.packages("h2o") }
  # # Next, we download packages that H2O depends on.
  # if (! ("methods" %in% rownames(installed.packages()))) { install.packages("methods") }
  # if (! ("statmod" %in% rownames(installed.packages()))) { install.packages("statmod") }
  # if (! ("stats" %in% rownames(installed.packages()))) { install.packages("stats") }
  # if (! ("graphics" %in% rownames(installed.packages()))) { install.packages("graphics") }
  # if (! ("RCurl" %in% rownames(installed.packages()))) { install.packages("RCurl") }
  # if (! ("jsonlite" %in% rownames(installed.packages()))) { install.packages("jsonlite") }
  # if (! ("tools" %in% rownames(installed.packages()))) { install.packages("tools") }
  # if (! ("utils" %in% rownames(installed.packages()))) { install.packages("utils") }
  # # Now we download, install and initialize the H2O package for R.
  # install.packages("h2o", type="source", repos=(c("http://h2o-release.s3.amazonaws.com/h2o/rel-tutte/2/R")))

  require("h2o")
  # set_all_stremr_options(estimator = "speedglm__glm")
  h2o::h2o.init(nthreads = -1)
  OData <- fitPropensity(OData,
                         gform_CENS = gform_CENS,
                         gform_TRT = gform_TRT,
                         models_TRT = defModel(estimator = "h2o__randomForest", ntrees = 20),
                         stratify_CENS = stratify_CENS)

})

