fit.h2oSuperLearner <- function(fit.class, fit, subsetH2Oframe, outvar, predvars, rows_subset, model_contrl, ...) {

  algorithms <- model_contrl$algorithm
  if (is.null(algorithms)) stop("must specify 'algorithm' name when running SuperLearner with 'h2o.grid'")
  if (!is.character(algorithms)) stop("'algorithm' must be a vector of strings naming the algorithms to use in 'h2o.grid'")

  fitted_models <- vector(mode = "list", length = length(algorithms))
  names(fitted_models) <- algorithms

  for (algorithm in algorithms) {
    grid_model_fit <- fit.h2ogrid(fit.class, fit, subsetH2Oframe, outvar, predvars, rows_subset, algorithm, model_contrl, ...)
    grid_model_H2O <- grid_model_fit$H2O.model.object
    fitted_models[[algorithm]] <- lapply(grid_model_H2O@model_ids, function(model_id) h2o::h2o.getModel(model_id))
  }

  # Put all of these models in a single list for stacking
  fitted_models_all <- NULL
  for (algorithm in algorithms) {
    fitted_models_all <- c(fitted_models_all, fitted_models[[algorithm]])
  }

  # Specify a defalt GLM as the metalearner
  metalearner <- model_contrl$metalearner
  stack <- h2oEnsemble::h2o.stack(models = fitted_models_all, response_frame = subsetH2Oframe[,outvar], metalearner = metalearner)

  # ----------------------------------------------------------------------------------------------------
  # Compute the final SL performance on the training set:
  # ----------------------------------------------------------------------------------------------------
  print("SuperLearner fit:"); print(stack$metafit)
  perf <- h2oEnsemble::h2o.ensemble_performance(stack, newdata = subsetH2Oframe, score_base_models = FALSE)
  print("SuperLearner overall performance (AUC) on the training set: "); print(perf)
  print("SuperLearner overall performance (MSE) on the training set: "); print(perf$ensemble@metrics$MSE)

  # h2o.glm_nn <- function(..., non_negative = TRUE) h2o.glm.wrapper(..., non_negative = non_negative)
  # stack3 <- h2o.metalearn(stack, metalearner = "h2o.glm_nn")
  # perf3 <- h2o.ensemble_performance(stack3, newdata = subsetH2Oframe, score_base_models = FALSE)
  # print(perf3)

  out_coef <- vector(mode = "numeric", length = length(fitted_models_all)+1)
  out_coef[] <- NA
  names(out_coef) <- names(stack$metafit@model$coefficients)
  out_coef[names(stack$metafit@model$coefficients)] <- stack$metafit@model$coefficients

  fit$coef <- out_coef;
  fit$linkfun <- NA
  fit$nobs <- length(rows_subset)

  if (gvars$verbose) {
    print("SuperLearner fits:")
    print(fit$coef)
  }

  fit$fitfunname <- "h2oEnsemble::h2o.stack";
  fit$H2O.model.object <- stack

  # browser()
  # TO DIRECTLY SAVE ALL MODELS FIT DURING GRID SEARCH
  # fit$fitted_models_all <- fitted_models_all
  class(fit) <- c(class(fit)[1], c("H2Oensemblemodel"))
  return(fit)
}