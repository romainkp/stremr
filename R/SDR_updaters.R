# SDR.updater.xgb <- function(x, y, offset, wts, params, ...) {
#   # params <- self$reg$SDR_model
#   if (gvars$verbose) { cat("running SDR w/ following params: \n "); str(params) }
#   xgb_dat <- xgboost::xgb.DMatrix(as.matrix(x), label = y)
#   xgboost::setinfo(xgb_dat, "base_margin", offset)
#   xgboost::setinfo(xgb_dat, "weight", wts)
#   nrounds <- params[["nrounds"]]
#   params[["nrounds"]] <- NULL

#   if (is.null(nrounds)) {
#     cat("...running cv to figure out best nrounds for epsilon target...\n")
#     mfitcv <- xgboost::xgb.cv(params = params, data = xgb_dat, nrounds = 100, nfold = 5, early_stopping_rounds = 10, verbose = 0)
#     nrounds <- mfitcv$best_iteration
#     cat("...best nrounds: ", nrounds, "\n")
#   }

#   if (gvars$verbose) { cat("...running SDR update with xgboost...\n") }
#   mfit <- xgboost::xgb.train(params = params, data = xgb_dat, nrounds = nrounds)

#   return(mfit)
# }

# SDR.updater.glm <- function(x, y, offset, wts, params, ...) {
#   cat("...running SDR update with GLM...\n")
#   mfit <- stats::glm.fit(x = cbind(Intercept = 1L, as.matrix(x)),
#                          y = y,
#                          weights = wts,
#                          offset = offset, # method=c('eigen','Cholesky','qr'),
#                          family = stats::quasibinomial(),
#                          control = glm.control(trace = FALSE))
#   return(mfit)
# }

# TMLE.updater.glm <- function(x, y, offset, wts, params, ...) {
#   cat("...running TMLE update with intercept-only GLM...\n")
#   mfit <- stats::glm.fit(x = matrix(1L, ncol = 1, nrow = length(y)),
#                          y = y,
#                          weights = wts,
#                          offset = offset, # method=c('eigen','Cholesky','qr'),
#                          family = stats::quasibinomial(),
#                          control = glm.control(trace = FALSE))
#   return(mfit)
# }

# SDR.predict.xgb <- function(mfit, x, offset, ...) {
#   xgb_dat <- xgboost::xgb.DMatrix(as.matrix(x))
#   xgboost::setinfo(xgb_dat, "base_margin", offset)
#   Qk_hat_star_all <- predict(mfit, xgb_dat)
#   return(Qk_hat_star_all)
# }

# SDR.predict.glm <- function(mfit, x, offset, ...) {
#   Xmat <- cbind(Intercept = 1L, as.matrix(x))
#   eta <- Xmat[,!is.na(mfit$coef), drop = FALSE] %*% mfit$coef[!is.na(mfit$coef)]
#   Qk_hat_star_all <- stremr:::logit_linkinv(offset + eta)
#   return(Qk_hat_star_all)
# }

# wrapper_to_SDRupdate_wrapper <- function(wrapper) {
#   if (is.character(wrapper)) {
#       print(wrapper)
#       wrapper <- eval(parse(text = wrapper))
#   }
#   function(Y, X, newX, ...) {
#       # nY <- ncol(Y)

#       offset <- X[, "offset"]
#       X <- X[, colnames(X)[!colnames(X) %in% "offset"], drop = FALSE]

#       new_offset <- newX[, "offset"]
#       newX <- newX[, colnames(newX)[!colnames(newX) %in% "offset"], drop = FALSE]

#       wrapper(Y, X, newX, offset, new_offset, ...)

#       col_fits <- lapply(seq_len(nY), function(col) {
#           wrapper(Y[, col], ...)
#       })

#       col_preds <- sapply(col_fits, `[[`, "pred")
#       fit <- list(col_fits = col_fits)
#       out <- list(pred = col_preds, fit = fit)


#       class(out$fit) <- c("mvSL.wrapper")
#       return(out)
#   }
# }

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
  fit.glm <- stats::glm.fit(x = matrix(1L, ncol = 1, nrow = length(Y)),
                         y = Y,
                         weights = obsWeights,
                         offset = offset, # method=c('eigen','Cholesky','qr'),
                         family = stats::quasibinomial(),
                         control = glm.control(trace = FALSE))

  fit <- list(object = list(est.fun = "stats::glm.fit", family = "quasibinomial", coef = fit.glm$coef))
  class(fit) <- "SDR.updater.TMLE"

  pred <- predict(fit, newX)
  out <- list(pred = pred, fit = fit)
  out
}

#' @export
SDR.updater.speedglmTMLE <- function(Y, X, newX, family, obsWeights, ...) {
  # cat("...running speedglm TMLE update with intercept-only GLM...\n")
  offset <- X[, "offset"]
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
SDR.updater.xgb <- function(Y, X, newX, family, obsWeights, params, ...) {
  # cat("...running SDR updater xgboost w/ following params: \n "); str(params)
  offset <- X[, "offset"]
  X <- X[, colnames(X)[!colnames(X) %in% "offset"], drop = FALSE]

  xgb_dat <- xgboost::xgb.DMatrix(as.matrix(X), label = Y)
  xgboost::setinfo(xgb_dat, "base_margin", offset)
  xgboost::setinfo(xgb_dat, "weight", obsWeights)
  nrounds <- params[["nrounds"]]
  params[["nrounds"]] <- NULL

  if (is.null(nrounds)) {
    # cat("...running cv to figure out best nrounds for epsilon target...\n")
    mfitcv <- xgboost::xgb.cv(params = params, data = xgb_dat, nrounds = 100, nfold = 5, early_stopping_rounds = 10, verbose = 0)
    nrounds <- mfitcv$best_iteration
    cat("...best nrounds: ", nrounds, "\n")
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
  pred <- stremr:::logit_linkinv(offset + mfit$coef)
  pred
}

#' @export
predict.SDR.updater.glm <- function(object, newdata, ...) {
  mfit <- object$object
  offset <- newdata[, "offset"]
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
  newdata <- newdata[, colnames(newdata)[!colnames(newdata) %in% "offset"], drop = FALSE]
  xgb_dat <- xgboost::xgb.DMatrix(as.matrix(newdata))
  xgboost::setinfo(xgb_dat, "base_margin", offset)
  pred <- predict(mfit, xgb_dat)
  pred
}





