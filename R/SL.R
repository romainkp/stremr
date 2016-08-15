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

  # metalearner <- "h2o.glm.wrapper"
  stack <- h2oEnsemble::h2o.stack(models = fitted_models_all, response_frame = subsetH2Oframe[,outvar], metalearner = metalearner)

  # str(stack)
  print("SuperLearner fit:"); print(stack$metafit)
  # browser()

  # Compute test set performance:
  perf <- h2oEnsemble::h2o.ensemble_performance(stack, newdata = subsetH2Oframe)
  print("SuperLearner overall performance on the training set: "); print(perf)
  # h2o.glm_nn <- function(..., non_negative = TRUE) h2o.glm.wrapper(..., non_negative = non_negative)
  # stack3 <- h2o.metalearn(stack, metalearner = "h2o.glm_nn")
  # perf3 <- h2o.ensemble_performance(stack3, newdata = subsetH2Oframe, score_base_models = FALSE)
  # print(perf3)

  fit$coef <- NULL;
  fit$fitfunname <- "h2o.SL";
  fit$H2O.model.object <- stack
  class(fit) <- c(class(fit)[1], c("H2Oensemblemodel"))
  return(fit)

    # Random Grid Search (e.g. 120 second maximum)
    # This is set to run fairly quickly, increase max_runtime_secs
    # or max_models to cover more of the hyperparameter space.
    # Also, you can expand the hyperparameter space of each of the
    # algorithms by modifying the hyper param code below.

    # search_criteria <- list(strategy = "RandomDiscrete", max_runtime_secs = 20)
    # search_criteria <- list(strategy = "RandomDiscrete", max_runtime_secs = 120)
    # search_criteria <- list(strategy = "RandomDiscrete", max_runtime_secs = 5*120)
    # search_criteria <- list(strategy = "RandomDiscrete", max_models = 42, max_runtime_secs = 28800)
    # search_criteria <- list(strategy = "RandomDiscrete", stopping_metric = "AUTO", stopping_tolerance = 0.001, stopping_rounds = 10)
    # search_criteria <- list(strategy = "RandomDiscrete", stopping_metric = "misclassification", stopping_tolerance = 0.00001, stopping_rounds = 5)
    # nfolds <- 5
    # nfolds <- 10

    # # ------------------------------------------------------------
    # # GBM Hyperparamters:
    # # ------------------------------------------------------------
    # ntrees_opt <- c(100, 200, 300, 500)
    # learn_rate_opt <- c(0.005, 0.01, 0.03, 0.06)
    # max_depth_opt <- c(3, 4, 5, 6, 9)
    # sample_rate_opt <- c(0.7, 0.8, 0.9, 1.0)
    # col_sample_rate_opt <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
    # balance_classes_opt <- c(TRUE, FALSE)

    # hyper_params <- list(ntrees = ntrees_opt,
    #                     learn_rate = learn_rate_opt,
    #                      max_depth = max_depth_opt,
    #                      sample_rate = sample_rate_opt,
    #                      col_sample_rate = col_sample_rate_opt,
    #                      balance_classes = balance_classes_opt)
    # gbm_grid <- h2o.grid("gbm", x = predvars, y = outvar,
    #                      training_frame = subsetH2Oframe,
    #                      seed = 1,
    #                      nfolds = nfolds,
    #                      fold_assignment = "Modulo",
    #                      keep_cross_validation_predictions = TRUE,
    #                      hyper_params = hyper_params,
    #                      search_criteria = search_criteria)
    # gbm_models <- lapply(gbm_grid@model_ids, function(model_id) h2o.getModel(model_id))

    # # ------------------------------------------------------------
    # # RF Hyperparamters
    # # ------------------------------------------------------------
    # ntrees_opt <- c(100, 200, 300, 500)
    # mtries_opt <- 8:20
    # max_depth_opt <- c(5, 10, 15, 20, 25)
    # sample_rate_opt <- c(0.7, 0.8, 0.9, 1.0)
    # col_sample_rate_per_tree_opt <- c(0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)
    # balance_classes_opt <- c(TRUE, FALSE)

    # hyper_params <- list(ntrees = ntrees_opt,
    #                      mtries = mtries_opt,
    #                      max_depth = max_depth_opt,
    #                      sample_rate = sample_rate_opt,
    #                      col_sample_rate_per_tree = col_sample_rate_per_tree_opt,
    #                      balance_classes = balance_classes_opt)

    # rf_grid <- h2o.grid("randomForest", x = predvars, y = outvar,
    #                     training_frame = subsetH2Oframe,
    #                     seed = 1,
    #                     nfolds = nfolds,
    #                     fold_assignment = "Modulo",
    #                     keep_cross_validation_predictions = TRUE,
    #                     hyper_params = hyper_params,
    #                     search_criteria = search_criteria)
    # rf_models <- lapply(rf_grid@model_ids, function(model_id) h2o.getModel(model_id))

    # # ------------------------------------------------------------
    # # Deeplearning Hyperparamters
    # # ------------------------------------------------------------
    # activation_opt <- c("Rectifier", "RectifierWithDropout",
    #                     "Maxout", "MaxoutWithDropout")
    # hidden_opt <- list(c(10,10), c(20,15), c(50,50,50))
    # l1_opt <- c(0, 1e-3, 1e-5)
    # l2_opt <- c(0, 1e-3, 1e-5)
    # hyper_params <- list(activation = activation_opt,
    #                      hidden = hidden_opt,
    #                      l1 = l1_opt,
    #                      l2 = l2_opt)

    # dl_grid <- h2o.grid("deeplearning", x = predvars, y = outvar,
    #                     training_frame = subsetH2Oframe,
    #                     epochs = 15,
    #                     seed = 1,
    #                     nfolds = nfolds,
    #                     fold_assignment = "Modulo",
    #                     keep_cross_validation_predictions = TRUE,
    #                     hyper_params = hyper_params,
    #                     search_criteria = search_criteria)
    # dl_models <- lapply(dl_grid@model_ids, function(model_id) h2o.getModel(model_id))

    # ------------------------------------------------------------
    # GLM Hyperparamters
    # ------------------------------------------------------------
    # alpha_opt <- c(0,1,seq(0.1,0.9,0.1))
    # lambda_opt <- c(0,1e-7,1e-5,1e-3,1e-1)
    # hyper_params <- list(alpha = alpha_opt, lambda = lambda_opt)
    # h2o.no_progress()
    # glm_grid <- h2o::h2o.grid("glm", x = predvars, y = outvar,
    #                      training_frame = subsetH2Oframe,
    #                      family = "binomial",
    #                      nfolds = nfolds,
    #                      fold_assignment = "Modulo",
    #                      keep_cross_validation_predictions = TRUE,
    #                      hyper_params = hyper_params,
    #                      search_criteria = search_criteria)
    # glm_models <- lapply(glm_grid@model_ids, function(model_id) h2o::h2o.getModel(model_id))

    # Create a list of all the base models
    # models <- c(gbm_models, rf_models, dl_models, glm_models)
    # models <- c(glm_models)
    # browser()
}