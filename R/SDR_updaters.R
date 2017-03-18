#' @export
SDR.updater.NULL <- function(Y, X, newX, family, obsWeights, ...) {
  # cat("...running no update (epsilon update set to 0)...\n")
  fit <- list(object = list(est.fun = "stats::glm.fit", family = "stats::quasibinomial", coef = 0))
  class(fit) <- "SDR.updater.TMLE"
  pred <- predict(fit, newX)
  out <- list(pred = pred, fit = fit)
  out
}

#' @export
SDR.updater.glmTMLE <- function(Y, X, newX, family, obsWeights, ...) {
  # cat("...running glm TMLE update with intercept-only GLM...\n")
  offset <- X[, "offset"]
  offset[offset == Inf] <- 20
  # offset[offset == -Inf] <- -20
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
  class(fit) <- "SDR.updater.TMLE"

  pred <- predict(fit, newX)
  out <- list(pred = pred, fit = fit)
  out
}

#' @export
SDR.updater.speedglmTMLE <- function(Y, X, newX, family, obsWeights, ...) {
  # cat("...running speedglm TMLE update with intercept-only GLM...\n")
  offset <- X[, "offset"]
  offset[offset == Inf] <- 20
  # offset[offset == -Inf] <- -20
  fit.glm <- speedglm::speedglm.wfit(X = matrix(1L, ncol = 1, nrow = length(Y)),
                                     y = Y, weights = obsWeights, offset = offset,
                                     # method=c('eigen','Cholesky','qr'),
                                     family = quasibinomial(), trace = FALSE, maxit = 1000)

  fit <- list(object = list(est.fun = "speedglm::speedglm.wfit", family = "quasibinomial", coef = fit.glm$coef))
  class(fit) <- "SDR.updater.TMLE"

  pred <- predict(fit, newX)
  out <- list(pred = pred, fit = fit)
  out
}

#' @export
SDR.updater.glm <- function(Y, X, newX, family, obsWeights, ...) {
  # cat("...running SDR update with GLM...\n")
  offset <- X[, "offset"]
  offset[offset == Inf] <- 20
  # offset[offset == -Inf] <- -20
  X <- X[, colnames(X)[!colnames(X) %in% "offset"], drop = FALSE]
  fit.glm <- stats::glm.fit(x = cbind(Intercept = 1L, as.matrix(X)),
                            y = Y,
                            weights = obsWeights,
                            offset = offset, # method=c('eigen','Cholesky','qr'),
                            family = stats::quasibinomial(),
                            control = glm.control(trace = FALSE))

  fit <- list(object = list(est.fun = "stats::glm.fit", family = "stats::quasibinomial", coef = fit.glm$coef))
  class(fit) <- "SDR.updater.glm"
  pred <- predict(fit, newX)
  out <- list(pred = pred, fit = fit)
  out
}

#' @export
SDR.updater.xgb.delta1 <- function(Y, X, newX, family, obsWeights, params, ...) {
  params[["max_delta_step"]] <- 1
  SDR.updater.xgb(Y, X, newX, family, obsWeights, params, ...)
}
#' @export
SDR.updater.xgb.delta2 <- function(Y, X, newX, family, obsWeights, params, ...) {
  params[["max_delta_step"]] <- 2
  SDR.updater.xgb(Y, X, newX, family, obsWeights, params, ...)
}
#' @export
SDR.updater.xgb.delta3 <- function(Y, X, newX, family, obsWeights, params, ...) {
  params[["max_delta_step"]] <- 3
  SDR.updater.xgb(Y, X, newX, family, obsWeights, params, ...)
}
#' @export
SDR.updater.xgb.delta4 <- function(Y, X, newX, family, obsWeights, params, ...) {
  params[["max_delta_step"]] <- 4
  SDR.updater.xgb(Y, X, newX, family, obsWeights, params, ...)
}

#' @export
SDR.updater.xgb <- function(Y, X, newX, family, obsWeights, params, ...) {
  # cat("...running SDR updater xgboost w/ following params: \n "); str(params)
  offset <- X[, "offset"]
  offset[offset == Inf] <- 20
  # offset[offset == -Inf] <- -20
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
  class(fit) <- "SDR.updater.xgb"

  pred <- predict(fit, newX, ...)
  out <- list(pred = pred, fit = fit)
  out
}

#' @export
predict.SDR.updater.TMLE <- function(object, newdata, ...) {
  mfit <- object$object
  offset <- newdata[, "offset"]
  offset[offset == Inf] <- 20
  # offset[offset == -Inf] <- -20
  pred <- stremr:::logit_linkinv(offset + mfit$coef)
  pred
}

#' @export
predict.SDR.updater.glm <- function(object, newdata, ...) {
  mfit <- object$object
  offset <- newdata[, "offset"]
  offset[offset == Inf] <- 20
  # offset[offset == -Inf] <- -20
  newdata <- newdata[, colnames(newdata)[!colnames(newdata) %in% "offset"], drop = FALSE]

  Xmat <- cbind(Intercept = 1L, as.matrix(newdata))
  eta <- Xmat[,!is.na(mfit$coef), drop = FALSE] %*% mfit$coef[!is.na(mfit$coef)]
  pred <- stremr:::logit_linkinv(offset + eta)
  pred
}

#' @export
predict.SDR.updater.xgb <- function (object, newdata, ...)  {
  mfit <- object$object
  offset <- newdata[, "offset"]
  offset[offset == Inf] <- 20
  # offset[offset == -Inf] <- -20
  newdata <- newdata[, colnames(newdata)[!colnames(newdata) %in% "offset"], drop = FALSE]
  xgb_dat <- xgboost::xgb.DMatrix(as.matrix(newdata))
  xgboost::setinfo(xgb_dat, "base_margin", offset)
  pred <- predict(mfit, xgb_dat)
  pred
}





