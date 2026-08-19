.nma_mbnma_dosefun <- function(model, degree = 2L, user_formula = NULL) {
  model <- match.arg(model, c("linear", "quadratic", "cubic", "emax", "exponential", "loglinear", "spline", "user"))
  switch(model,
    linear = MBNMAdose::dpoly(degree = 1),
    quadratic = MBNMAdose::dpoly(degree = 2),
    cubic = MBNMAdose::dpoly(degree = 3),
    emax = MBNMAdose::demax(),
    exponential = MBNMAdose::dexp(),
    loglinear = MBNMAdose::dloglin(),
    spline = MBNMAdose::dspline(),
    user = {
      if (is.null(user_formula)) stop("`user_formula` is required for model='user'.", call. = FALSE)
      MBNMAdose::duser(fun = user_formula)
    }
  )
}

#' Fit dose-response network meta-analysis
#'
#' @param data For Bayesian MBNMA, an arm-level data frame in MBNMAdose format.
#'   For frequentist analysis, a contrast-level data frame accepted by `netdose::netdose()`.
#' @param framework `bayesian` or `frequentist`.
#' @param model Dose-response model. Frequentist choices include linear, exponential,
#'   quadratic, restricted cubic spline and fractional polynomials. Bayesian choices
#'   include polynomial, Emax, exponential, log-linear, spline and user-defined forms.
#' @param method Common/random model setting for MBNMAdose; for netdose use `random`/`common` in `...`.
#' @param user_formula Optional MBNMAdose `duser()` formula.
#' @param ... Passed to `mbnma.network()`/`mbnma.run()` or `netdose::netdose()`.
#' @return An `nmaflow_fit` object.
#' @export
nma_dose <- function(data, framework = c("bayesian", "frequentist"),
                     model = "emax", method = "random", user_formula = NULL, ...) {
  framework <- match.arg(framework)
  if (framework == "frequentist") {
    .nma_require("netdose", "frequentist dose-response NMA")
    args <- c(list(data = data, method = model), list(...))
    fit <- do.call(netdose::netdose, args)
    return(.nma_new("nmaflow_fit", engine = "netdose", framework = "frequentist_dose",
                    fit = fit, data = data, specification = list(model = model)))
  }
  .nma_require("MBNMAdose", "Bayesian dose-response NMA")
  net <- MBNMAdose::mbnma.network(data)
  fun <- .nma_mbnma_dosefun(model, user_formula = user_formula)
  fit <- MBNMAdose::mbnma.run(net, fun = fun, method = method, ...)
  .nma_new("nmaflow_fit", engine = "MBNMAdose", framework = "bayesian_dose",
           fit = fit, network = net, data = data, specification = list(model = model, method = method))
}

#' Create or fit agronomically useful nonlinear dose-response specifications
#'
#' @param data Arm-level MBNMAdose data.
#' @param model `mitscherlich`, `quadratic`, `emax`, or `custom`.
#' @param run If `FALSE`, return only the model specification. If `TRUE`, fit with MBNMAdose.
#' @param custom_formula Custom MBNMAdose dose-response formula when `model='custom'`.
#' @param ... Passed to [nma_dose()].
#' @return Specification list or fitted model.
#' @export
nma_nonlinear <- function(data, model = c("mitscherlich", "quadratic", "emax", "custom"),
                          run = TRUE, custom_formula = NULL, ...) {
  model <- match.arg(model)
  formula <- switch(model,
    mitscherlich = ~ beta.1 * (1 - exp(-beta.2 * dose)),
    quadratic = NULL,
    emax = NULL,
    custom = custom_formula
  )
  spec <- list(model = model, backend = "MBNMAdose",
               formula = formula,
               note = "Nonlinear parameters must be checked for identifiability and biological plausibility within the observed dose range.")
  if (!isTRUE(run)) return(spec)
  if (model == "quadratic") return(nma_dose(data, framework = "bayesian", model = "quadratic", ...))
  if (model == "emax") return(nma_dose(data, framework = "bayesian", model = "emax", ...))
  if (is.null(formula)) stop("A custom formula is required.", call. = FALSE)
  nma_dose(data, framework = "bayesian", model = "user", user_formula = formula, ...)
}

#' Predict dose-response curves
#' @param fit Dose-response `nmaflow_fit`.
#' @param ... Passed to the backend prediction method.
#' @export
nma_dose_curve <- function(fit, ...) {
  if (!inherits(fit, "nmaflow_fit")) stop("Expected nmaflow_fit.", call. = FALSE)
  if (fit$engine == "MBNMAdose") return(stats::predict(fit$fit, ...))
  if (fit$engine == "netdose") return(stats::predict(fit$fit, ...))
  stop("Dose prediction is not available for this backend.", call. = FALSE)
}

#' Estimate a finite-grid agronomic optimum from dose predictions
#'
#' @param fit Dose-response fit.
#' @param doses Numeric dose grid.
#' @param prediction_function Optional function returning a numeric response at each dose.
#' @param objective `max` or `min`.
#' @param ... Passed to prediction when a custom function is not used.
#' @return Dose and predicted response at the best finite-grid point.
#' @export
nma_dose_optimum <- function(fit, doses, prediction_function = NULL,
                             objective = c("max", "min"), ...) {
  objective <- match.arg(objective)
  if (!is.numeric(doses) || length(doses) < 2L) stop("Provide a numeric dose grid.", call. = FALSE)
  if (is.null(prediction_function)) {
    stop("Supply `prediction_function(fit, doses, ...)` because backend prediction structures differ by model.", call. = FALSE)
  }
  y <- prediction_function(fit, doses, ...)
  if (length(y) != length(doses)) stop("Prediction function must return one value per dose.", call. = FALSE)
  i <- if (objective == "max") which.max(y) else which.min(y)
  list(dose = doses[i], predicted = y[i], objective = objective,
       boundary = i %in% c(1L, length(doses)),
       note = "Finite-grid optimum; do not extrapolate beyond the experimental dose range without justification.")
}

#' Compare dose-response candidate models
#' @param fits Named list of dose-response fits.
#' @return Data frame of available information criteria.
#' @export
nma_dose_compare <- function(fits) {
  nma_regression_compare(fits)
}
