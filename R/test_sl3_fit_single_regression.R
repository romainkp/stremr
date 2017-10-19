
test_sl3_fit_single_regression <- function(data, nodes, models, predvars, outvar, subset_idx, ...) {

    folds <- data$make_origami_fold_from_column(subset_idx)
    task <- sl3::sl3_Task$new(data$dat.sVar[subset_idx, ],
                              covariates = predvars,
                              outcome = outvar,
                              id = data$nodes$IDnode,
                              folds = folds)

    reg_t <- system.time(
      model.fit <- try({models$train(task)})
      )
    print("reg_t"); print(reg_t)

    if (inherits(model.fit, "try-error")) {
      cat("\nsl3 error debugging info:\n");
      print(model.fit)
    }

  return(model.fit)
}