
truncate_offset <- function(offset,
                            up_trunc_offset = stremrOptions("up_trunc_offset"),
                            low_trunc_offset = stremrOptions("low_trunc_offset")) {
  offset[offset >= up_trunc_offset] <- up_trunc_offset
  offset[offset <= low_trunc_offset] <- low_trunc_offset
  return(offset)
}

#' iTMLE NULL learner / updater
#'
#' Returns the inputs Y as they are, without updating. This function is passed along as a separate learner
#' to the SuperLearner implementation of origami package
#' @param Y Input outcomes
#' @param X Input design matrix with training data.
#' Must contain a column named "offset", which contains the offsets converted to logit-linear scale.
#' @param newX Input design matrix with test data.
#' Same requirement as for \code{X}: must contain a column named "offset",
#' which contains the offsets converted to logit-linear scale.
#' @param family Link function (ignored).
#' @param obsWeights Row-specific weights
#' @param ... Additional arguments to be passed on to \code{origami} package.
#' @export
TMLE.updater.NULL <- function(Y, X, newX, family, obsWeights, ...) {
  # cat("...running no update (epsilon update set to 0)...\n")
  fit <- list(object = list(est.fun = "stats::glm.fit", family = "stats::quasibinomial", coef = 0))
  class(fit) <- "TMLE.updater"
  pred <- predict(fit, newX)
  out <- list(pred = pred, fit = fit)
  out
}

#' iTMLE univariate glm learner / updater
#'
#' Performs a univariate TMLE update. This function is passed along as a separate learner
#' to the SuperLearner implementation of origami package
#' @param Y Input outcomes
#' @param X Input design matrix with training data.
#' Must contain a column named "offset", which contains the offsets converted to logit-linear scale.
#' @param newX Input design matrix with test data.
#' Same requirement as for \code{X}: must contain a column named "offset",
#' which contains the offsets converted to logit-linear scale.
#' @param family Link function (ignored).
#' @param obsWeights Row-specific weights
#' @param ... Additional arguments to be passed on to \code{origami} package.
#' @export
TMLE.updater.glm <- function(Y, X, newX, family, obsWeights, ...) {
  eps_tol = stremrOptions("eps_tol")
  if (any(Y < 0) | any(Y > 1)) stop("TMLE update with glm cannot be performed for outcomes outside of 0-1.")
  Y[Y < eps_tol] <- eps_tol
  Y[Y > (1 - eps_tol)] <- (1 - eps_tol)

  # cat("...running glm TMLE update with intercept-only GLM...\n")
  offset <- truncate_offset(X[, "offset"])
  fit.glm <- try(stats::glm.fit(x = matrix(1L, ncol = 1, nrow = length(Y)),
                                 y = Y,
                                 weights = obsWeights,
                                 offset = offset, # method=c('eigen','Cholesky','qr'),
                                 family = stats::quasibinomial(),
                                 control = glm.control(trace = FALSE)))

  if (inherits(fit.glm, "try-error")) { # TMLE update w/ glm failed
    message("GLM TMLE update has failed during SDR, setting epsilon update to 0")
    warning("GLM TMLE update has failed during SDR, setting epsilon update to 0")
    epsilon_coef <- 0.0
  } else {
    epsilon_coef <- fit.glm$coef
  }

  fit <- list(object = list(est.fun = "stats::glm.fit", family = "quasibinomial", coef = epsilon_coef))
  class(fit) <- "TMLE.updater"

  pred <- predict(fit, newX)
  out <- list(pred = pred, fit = fit)
  out
}

#' TMLE glm learner / updater
#'
#' Performs TMLE update using a main-terms GLM (logistic regression with \code{speedglm} package).
#' This function may be passed to fit_GCOMP / fit_TMLE function.
#' May be also used as a separate learner for iTMLE.
#' @rdname TMLE.updater.glm
#' @export
TMLE.updater.speedglm <- function(Y, X, newX, family, obsWeights, ...) {
  eps_tol = stremrOptions("eps_tol")
  if (any(Y < 0) | any(Y > 1)) stop("TMLE update with speedglm cannot be performed for outcomes outside of 0-1.")
  Y[Y < eps_tol] <- eps_tol
  Y[Y > (1 - eps_tol)] <- (1 - eps_tol)

  # cat("...running speedglm TMLE update with intercept-only GLM...\n")
  offset <- truncate_offset(X[, "offset"])
  fit.glm <- speedglm::speedglm.wfit(X = matrix(1L, ncol = 1, nrow = length(Y)),
                                     y = Y, weights = obsWeights, offset = offset,
                                     # method=c('eigen','Cholesky','qr'),
                                     family = quasibinomial(), trace = FALSE, maxit = 1000)

  fit <- list(object = list(est.fun = "speedglm::speedglm.wfit", family = "quasibinomial", coef = fit.glm$coef))
  class(fit) <- "TMLE.updater"

  pred <- predict(fit, newX)
  out <- list(pred = pred, fit = fit)
  out
}

#' Linear TMLE glm updater
#'
#' Performs linear TMLE update using main-terms GLM (family=gaussian()), based on \code{speedglm} package.
#' This function may be passed to fit_GCOMP / fit_TMLE function.
#' May be also used as a separate learner for iTMLE.
#' @rdname TMLE.updater.glm
#' @export
linear.TMLE.updater.speedglm <- function(Y, X, newX, family, obsWeights, ...) {
  eps_tol = stremrOptions("eps_tol")
  if (any(Y < 0) | any(Y > 1)) stop("TMLE update with speedglm cannot be performed for outcomes outside of 0-1.")
  Y[Y < eps_tol] <- eps_tol
  Y[Y > (1 - eps_tol)] <- (1 - eps_tol)

  # cat("...running speedglm TMLE update with intercept-only GLM...\n")
  offset <- truncate_offset(plogis(X[, "offset"]))
  fit.glm <- speedglm::speedglm.wfit(X = matrix(1L, ncol = 1, nrow = length(Y)),
                                     y = Y, weights = obsWeights, offset = offset,
                                     # method=c('eigen','Cholesky','qr'),
                                     family = gaussian(), trace = FALSE, maxit = 1000)

  fit <- list(object = list(est.fun = "speedglm::speedglm.wfit", family = "gaussian", coef = fit.glm$coef))
  class(fit) <- "linear.TMLE.updater"

  pred <- predict(fit, newX)
  out <- list(pred = pred, fit = fit)
  out
}

#' iTMLE glm learner / updater
#'
#' Performs an SDR update using a main-terms GLM (logistic regression with \code{stats::glm}).
#' This function is passed along as a separate learner
#' to the SuperLearner implementation in \code{origami} package
#' @param Y Input outcomes
#' @param X Input design matrix with training data.
#' Must contain a column named "offset", which contains the offsets converted to logit-linear scale.
#' @param newX Input design matrix with test data.
#' Same requirement as for \code{X}: must contain a column named "offset",
#' which contains the offsets converted to logit-linear scale.
#' @param family Link function (ignored).
#' @param obsWeights Row-specific weights
#' @param ... Additional arguments to be passed on to \code{origami} package.
#' @export
iTMLE.updater.glm <- function(Y, X, newX, family, obsWeights, ...) {
  eps_tol = stremrOptions("eps_tol")
  if (any(Y < 0) | any(Y > 1)) stop("TMLE update with glm cannot be performed for outcomes outside of 0-1.")
  Y[Y < eps_tol] <- eps_tol
  Y[Y > (1 - eps_tol)] <- (1 - eps_tol)

  # cat("...running SDR update with GLM...\n")
  offset <- truncate_offset(X[, "offset"])
  X <- X[, colnames(X)[!colnames(X) %in% "offset"], drop = FALSE]
  fit.glm <- stats::glm.fit(x = cbind(Intercept = 1L, as.matrix(X)),
                            y = Y,
                            weights = obsWeights,
                            offset = offset, # method=c('eigen','Cholesky','qr'),
                            family = stats::quasibinomial(),
                            control = glm.control(trace = FALSE))

  fit <- list(object = list(est.fun = "stats::glm.fit", family = "stats::quasibinomial", coef = fit.glm$coef))
  class(fit) <- "iTMLE.updater.glm"
  pred <- predict(fit, newX)
  out <- list(pred = pred, fit = fit)
  out
}

#' iTMLE gbm learner / updater
#'
#' Performs an SDR update using GBM from \code{xgboost} package with parameter \code{max_delta_step}=1.
#' This function is passed along as a separate learner
#' to the SuperLearner implementation of origami package
#' @rdname iTMLE.updater.xgb
#' @export
iTMLE.updater.xgb.delta1 <- function(Y, X, newX, family, obsWeights, params, ...) {
  params[["max_delta_step"]] <- 1
  iTMLE.updater.xgb(Y, X, newX, family, obsWeights, params, ...)
}

#' iTMLE gbm learner / updater
#'
#' Performs an SDR update using GBM from \code{xgboost} package with parameter \code{max_delta_step}=2.
#' This function is passed along as a separate learner
#' to the SuperLearner implementation of origami package
#' @rdname iTMLE.updater.xgb
#' @export
iTMLE.updater.xgb.delta2 <- function(Y, X, newX, family, obsWeights, params, ...) {
  params[["max_delta_step"]] <- 2
  iTMLE.updater.xgb(Y, X, newX, family, obsWeights, params, ...)
}

#' iTMLE gbm learner / updater
#'
#' Performs an SDR update using GBM from \code{xgboost} package with parameter \code{max_delta_step}=3.
#' This function is passed along as a separate learner
#' to the SuperLearner implementation of origami package
#' @rdname iTMLE.updater.xgb
#' @export
iTMLE.updater.xgb.delta3 <- function(Y, X, newX, family, obsWeights, params, ...) {
  params[["max_delta_step"]] <- 3
  iTMLE.updater.xgb(Y, X, newX, family, obsWeights, params, ...)
}

#' iTMLE gbm learner / updater
#'
#' Performs an SDR update using GBM from \code{xgboost} package with parameter \code{max_delta_step}=4.
#' This function is passed along as a separate learner
#' to the SuperLearner implementation of origami package
#' @rdname iTMLE.updater.xgb
#' @export
iTMLE.updater.xgb.delta4 <- function(Y, X, newX, family, obsWeights, params, ...) {
  params[["max_delta_step"]] <- 4
  iTMLE.updater.xgb(Y, X, newX, family, obsWeights, params, ...)
}

#' iTMLE gbm learner / updater
#'
#' General GBM learner from \code{xgboost} package.
#' This function is passed along as a separate learner
#' to the SuperLearner implementation of origami package
#' @param Y Input outcomes
#' @param X Input design matrix with training data.
#' Must contain a column named "offset", which contains the offsets converted to logit-linear scale.
#' @param newX Input design matrix with test data.
#' Same requirement as for \code{X}: must contain a column named "offset",
#' which contains the offsets converted to logit-linear scale.
#' @param family Link function (ignored).
#' @param obsWeights Row-specific weights
#' @param params Tuning parameters passed on to \code{xgboost}.
#' @param ... Additional arguments to be passed on to \code{origami} package.
#' @export
iTMLE.updater.xgb <- function(Y, X, newX, family, obsWeights, params, ...) {
  # cat("...running SDR updater xgboost w/ following params: \n "); str(params)
  offset <- truncate_offset(X[, "offset"])
  X <- X[, colnames(X)[!colnames(X) %in% "offset"], drop = FALSE]

  xgb_dat <- xgboost::xgb.DMatrix(as.matrix(X), label = Y)
  xgboost::setinfo(xgb_dat, "base_margin", offset)
  xgboost::setinfo(xgb_dat, "weight", obsWeights)
  nrounds <- params[["nrounds"]]
  params[["nrounds"]] <- NULL

  if (is.null(nrounds)) {
    cat("...running cv to figure out best nrounds for epsilon target...\n")
    mfitcv <- xgboost::xgb.cv(params = params, data = xgb_dat, nrounds = 100, nfold = 5, early_stopping_rounds = 10, verbose = 0)
    nrounds <- mfitcv$best_iteration
    # cat("...best nrounds: ", nrounds, "\n")
  }

  fit.xgb <- xgboost::xgb.train(params = params, data = xgb_dat, nrounds = nrounds)
  fit <- list(object = fit.xgb)
  class(fit) <- "iTMLE.updater.xgb"

  pred <- predict(fit, newX, ...)
  out <- list(pred = pred, fit = fit)
  out
}

#' @param object Results of calling \code{TMLE.updater.speedglm} or \code{TMLE.updater.glm}.
#' @param newdata Design matrix with test data for which predictions should be obtained.
#' Must contain a column named "offset".
#' @rdname TMLE.updater.glm
#' @export
predict.TMLE.updater <- function(object, newdata, offset, ...) {
  mfit <- object$object
  if (missing(offset)) {
    offset <- newdata[, "offset"]
  }
  offset <- truncate_offset(offset)
  pred <- logit_linkinv(offset + mfit$coef)
  pred
}

#' @param object Results of calling \code{linear.TMLE.updater.speedglm}.
#' @param newdata Design matrix with test data for which predictions should be obtained.
#' Must contain a column named "offset".
#' @rdname linear.TMLE.updater.speedglm
#' @export
predict.linear.TMLE.updater <- function(object, newdata, offset, ...) {
  mfit <- object$object
  if (missing(offset)) {
    offset <- newdata[, "offset"]
  }
  offset <- truncate_offset(offset)
  pred <- offset + mfit$coef
  pred
}

#' @param object Results of calling \code{iTMLE.updater.glm}.
#' @param newdata Design matrix with test data for which predictions should be obtained.
#' Must contain a column named "offset".
#' @rdname iTMLE.updater.glm
#' @export
predict.iTMLE.updater.glm <- function(object, newdata, offset, ...) {
  mfit <- object$object
  if (missing(offset)) {
    offset <- newdata[, "offset"]
  }
  offset <- truncate_offset(offset)
  newdata <- newdata[, colnames(newdata)[!colnames(newdata) %in% "offset"], drop = FALSE]
  Xmat <- cbind(Intercept = 1L, as.matrix(newdata))
  eta <- Xmat[,!is.na(mfit$coef), drop = FALSE] %*% mfit$coef[!is.na(mfit$coef)]
  pred <- logit_linkinv(offset + eta)
  pred
}

#' @param object Results of calling \code{iTMLE.updater.xgb} functions.
#' @param newdata Design matrix with test data for which predictions should be obtained.
#' Must contain a column named "offset".
#' @rdname iTMLE.updater.xgb
#' @export
predict.iTMLE.updater.xgb <- function (object, newdata, offset, ...)  {
  mfit <- object$object
  if (missing(offset)) {
    offset <- newdata[, "offset"]
  }
  offset <- truncate_offset(offset)
  newdata <- newdata[, colnames(newdata)[!colnames(newdata) %in% "offset"], drop = FALSE]
  xgb_dat <- xgboost::xgb.DMatrix(as.matrix(newdata))
  xgboost::setinfo(xgb_dat, "base_margin", offset)
  pred <- predict(mfit, xgb_dat)
  pred
}





