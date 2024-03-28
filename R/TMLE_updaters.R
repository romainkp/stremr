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

  if (sum(abs(obsWeights)) < 10^-9) {
    message("GLM TMLE update cannot be performed since all IP-weights are exactly zero, setting epsilon to 0")
    epsilon <- 0
  } else if (all(Y[obsWeights > 0] < eps_tol) || all(Y[obsWeights > 0] > (1 - eps_tol))) {
    message("GLM TMLE update cannot be performed since the outcomes (Y) are either all 0 or all 1, setting epsilon to 0")
    epsilon <- 0
  } else {
    Y[Y < eps_tol] <- eps_tol
    Y[Y > (1 - eps_tol)] <- (1 - eps_tol)

    # cat("...running speedglm TMLE update with intercept-only GLM...\n")
    offset <- truncate_offset(X[, "offset"])
    fit.glm <- try(speedglm::speedglm.wfit(X = matrix(1L, ncol = 1, nrow = length(Y)),
                                         y = Y, weights = obsWeights, offset = offset,
                                         # method=c('eigen','Cholesky','qr'),
                                         family = quasibinomial(), trace = FALSE, maxit = 1000),
                silent = TRUE)

    if (inherits(fit.glm, "try-error")) { # TMLE update failed
      # if (gvars$verbose)
      message("GLM TMLE update has failed, setting epsilon = 0")
      # warning("GLM TMLE update has failed, setting epsilon = 0")
      epsilon <- 0
    } else {
      epsilon <- fit.glm$coef
    }
  }

  if (gvars$verbose) { cat("\ncoef: ", epsilon,"\n") }

  fit <- list(object = list(est.fun = "speedglm::speedglm.wfit", family = "quasibinomial", coef = epsilon))
  class(fit) <- "TMLE.updater"

  pred <- predict(fit, newX)
  out <- list(pred = pred, fit = fit)
  out
}

#' TMLE glm learner / updater for rare outcomes
#'
#' Performs TMLE update using a main-terms logistic regression when the outcome is not between 0 and 1 but yet its conditional mean is between 0 and 1.
#' In that case glm cannot be used to fit the logistic regression as it will return an error but because the negative quasi-log-likelihood remains a valid loss function for the outcome then a logistic regression fit can be optained by minimizing the loss manually using an optimization function like optim.
#' This rare outcome updater function may be passed to fit_GCOMP / fit_TMLE function.
#' May be also used as a separate learner for iTMLE.
#' This implements the rare outcomes functions by Laura Balzar shared here: https://www.dropbox.com/s/bs38e8mjt549kgo/rTMLE_sample_code.R?dl=0 (also stored on the ON TARGET DM gitlab repository)
#' specifically: rtmle.wt(), tmleLogLike(), tmlegrad() 
#' Also, see "Estimating Effects with Rare Outcomes and High Dimensional Covariates: Knowledge is Power" paper here: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5436729/
#' @rdname TMLE.updater.glm
#' @export
TMLE.updater.rare <- function(Y, X, newX, family, obsWeights, ...) {

  ## extract stremrOptions values
  eps_tol = stremrOptions("eps_tol") # default 1e-05
  up_rare_y = stremrOptions("up_rare_y") # default 0.05
  low_rare_y = stremrOptions("low_rare_y") # default 0.00
  
  ## checks
  if (any(Y < 0) | any(Y > 1)) stop("\nRare Outcome TMLE update cannot be performed for outcomes outside of 0-1.")
  
  ## rare outcomes implemented in last else{}
  if (sum(abs(obsWeights)) < 10^-9) {
    message("\nRare Outcome TMLE update cannot be performed since all IP-weights are exactly zero, setting epsilon to 0")
    epsilon <- 0
    # define fit / pred when weights are 0
    fit <- list(object = list(est.fun = "speedglm::speedglm.wfit", family = "quasibinomial", coef = epsilon))
    class(fit) <- "TMLE.updater"
    pred <- predict(fit, newX)
  } else if (all(Y[obsWeights > 0] < eps_tol) || all(Y[obsWeights > 0] > (1 - eps_tol))) {
    message("\nRare Outcome TMLE update cannot be performed since the outcomes (Y) are either all 0 or all 1, setting epsilon to 0")
    epsilon <- 0
    # define fit / pred when Y is constant
    fit <- list(object = list(est.fun = "speedglm::speedglm.wfit", family = "quasibinomial", coef = epsilon))
    class(fit) <- "TMLE.updater"
    pred <- predict(fit, newX)
  } else {
    # update Y values to be between 0 and 1 instead of 0 and 1 
    Y[Y < eps_tol] <- eps_tol
    Y[Y > (1 - eps_tol)] <- (1 - eps_tol)
    
    ## tmleLogLike() - updated from Laura's code, removed clever covariates
    tmleLogLike<- function(eps, Y, Q.AW, wt){
      pi <- plogis( qlogis(Q.AW) + eps)
      pi[pi==0] <- .Machine$double.neg.eps		
      pi[pi==1] <- 1-.Machine$double.neg.eps
      logLike <- sum( wt*( Y*log(pi)+ (1-Y)*log(1-pi)) )
      
      if(is.na(logLike) | is.nan(logLike) | is.infinite(logLike)){
        eps <- 0	
        logLike <- 99
      }
      return(-logLike)		
    }
    ## tmlegrad() - updated from Laura's code, removed clever covariates   
    tmlegrad <- function(eps, Y, Q.AW, wt){
      
      pi <- plogis( qlogis(Q.AW) + eps)
      pi[pi==0] <- .Machine$double.neg.eps
      pi[pi==1]<- 1-.Machine$double.neg.eps
      resid <- wt*(Y-pi) # weighted residuals
      gr <- crossprod( cbind(rep(1, length(resid))), resid)
      
      if( sum(is.na(gr))>0 | sum(is.nan(gr))>0 | sum(is.infinite(gr))>0 ){
        gr <- 99
      }
      return(-gr)		
    }
    
    ## rtmle.wt() - steps / implementation
    ## STEP 0 of  - transform outcome --------------------------------------
    Y.tilde <- (Y-low_rare_y)/(up_rare_y-low_rare_y)
    
    ## STEP 1 - initial estimate --------------------------------------
    ## NOTE: we fit the initial estimate with fit_GCOMP and glmnet with the original Y (not transformed Y.tilde as Laura does using optim())
    ## Because Laura fit the optim() in step 1 w/ transformed Y.tilde, we must transform the offset to match the scale of the offset Laura 
    ## gets from her step 1 before moving to step 2. Laura's offset is Q.opt$QbarAW.
    
    ### update low_trunc_offset and up_trunc_offset
    # stremrOptions("up_trunc_offset", -log((1-min((up_rare_y - eps_tol), up_rare_y))/(min((up_rare_y - eps_tol), up_rare_y))))
    stremrOptions("up_trunc_offset", -log((1-up_rare_y)/(up_rare_y)))
    stremrOptions("low_trunc_offset", -log((1-max(low_rare_y,.Machine$double.neg.eps))/(max(low_rare_y,.Machine$double.neg.eps))))

    log.offset <- truncate_offset(X[, "offset"]) ## truncate offset, on log scale, range [low_trunc_offset,up_trunc_offset]
    prob.offset <- 1/(1+exp(-log.offset)) ## transform to probability scale, range [low_rare_y, up_rare_y]
    offset.tilde <- (prob.offset-low_rare_y)/(up_rare_y-low_rare_y) ## used in optim() in step 2, 
                                                                    ## same scale as Laura's Q.opt$QbarAW
                                                                    ## range [0.0, 1.0]
    log.offset.tilde <- log(offset.tilde/(1-offset.tilde)) ## used in predict() in step 3, will be NA if offset.tilde is 1.0, 
    log.offset.tilde[is.na(log.offset.tilde)] <- Inf       ## Replace NA values with Inf bc occurs when #/(1-1) which is Inf                                                      


    ## STEP 2 - targeting --------------------------------------
    ## NOTE: Did not implement parametric estimate / clever covariates code from Laura
    ## Already ran our PS and have our weights
    
    # initial parameter estimate is 0 (par)
    fit.glm <- optim(par=0, fn=tmleLogLike, gr=tmlegrad, Y=Y.tilde, Q.AW=offset.tilde, 
                     wt=obsWeights,  method="BFGS", hessian=F,  control = list(maxit = 500, trace=F)) 
    
    if (inherits(fit.glm, "try-error")) { # TMLE update failed
      if (gvars$verbose)
        message("\nRare Outcome TMLE update has failed, setting epsilon = 0")
      warning("\nGLM TMLE update has failed, setting epsilon = 0")
      epsilon <- 0
    } else {
      # assign epsilon
      epsilon <- fit.glm$par
    }
    
    # fit rare outcomes TMLE and output pred with transformation
    fit <- list(object = list(est.fun = "rare.outcome.w.optim", family = "quasibinomial", coef = epsilon))
    class(fit) <- "TMLE.updater"
    pred <- predict(fit, newX, log.offset.tilde)
    pred <- low_rare_y + (up_rare_y-low_rare_y)*pred # range [low_rare_y, up_rare_y]
    
    if (gvars$verbose) {
      message("\nRare Outcome Pred: printing summary(pred) after transformation")
      print(summary(pred))
    }
  }
  
  if (gvars$verbose) { cat("\ncoef: ", epsilon,"\n") }
  
  out <- list(pred = pred, fit = fit) 
  out
}


#' Linear TMLE glm updater
#'
#' Performs linear TMLE update using main-terms GLM (family=gaussian()), based on \code{speedglm} package.
#' This function may be passed to fit_GCOMP / fit_TMLE function.
#' May be also used as a separate learner for iTMLE.
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
#' May contain a column named "offset", instead of it being passed as a separate argument.
#' @param offset Offset (on logit scale if using logistic regression update).
#' @rdname TMLE.updater.glm
#' @export
predict.TMLE.updater <- function(object, newdata, offset, ...) {
  mfit <- object$object
  if (missing(offset)) {
    offset <- newdata[, "offset"]
  }
  # offset <- truncate_offset(offset)
  pred <- plogis(offset + mfit$coef)
  pred
}

#' @param object Results of calling \code{linear.TMLE.updater.speedglm}.
#' @param newdata Design matrix with test data for which predictions should be obtained.
#' May contain a column named "offset", instead of it being passed as a separate argument.
#' @param offset Offset (on logit scale if using logistic regression update).
#' @rdname linear.TMLE.updater.speedglm
#' @export
predict.linear.TMLE.updater <- function(object, newdata, offset, ...) {
  mfit <- object$object
  if (missing(offset)) {
    offset <- newdata[, "offset"]
  }
  # offset <- truncate_offset(offset)
  pred <- offset + mfit$coef
  pred
}

#' @param object Results of calling \code{iTMLE.updater.glm}.
#' @param newdata Design matrix with test data for which predictions should be obtained.
#' May contain a column named "offset", instead of it being passed as a separate argument.
#' @param offset Offset (on logit scale if using logistic regression update).
#' @rdname iTMLE.updater.glm
#' @export
predict.iTMLE.updater.glm <- function(object, newdata, offset, ...) {
  mfit <- object$object
  if (missing(offset)) {
    offset <- newdata[, "offset"]
  }
  # offset <- truncate_offset(offset)
  newdata <- newdata[, colnames(newdata)[!colnames(newdata) %in% "offset"], drop = FALSE]
  Xmat <- cbind(Intercept = 1L, as.matrix(newdata))
  eta <- Xmat[,!is.na(mfit$coef), drop = FALSE] %*% mfit$coef[!is.na(mfit$coef)]
  pred <- plogis(offset + eta)
  pred
}

#' @param object Results of calling \code{iTMLE.updater.xgb} functions.
#' @param newdata Design matrix with test data for which predictions should be obtained.
#' May contain a column named "offset", instead of it being passed as a separate argument.
#' @param offset Offset (on logit scale if using logistic regression update).
#' @rdname iTMLE.updater.xgb
#' @export
predict.iTMLE.updater.xgb <- function (object, newdata, offset, ...)  {
  mfit <- object$object
  if (missing(offset)) {
    offset <- newdata[, "offset"]
  }
  # offset <- truncate_offset(offset)
  newdata <- newdata[, colnames(newdata)[!colnames(newdata) %in% "offset"], drop = FALSE]
  xgb_dat <- xgboost::xgb.DMatrix(as.matrix(newdata))
  xgboost::setinfo(xgb_dat, "base_margin", offset)
  pred <- predict(mfit, xgb_dat)
  pred
}





