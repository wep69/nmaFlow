#' Fit component network meta-analysis
#'
#' @param fit A standard netmeta-backed `nmaflow_fit` where complex treatment names encode
#'   components using a separator.
#' @param inactive Optional inactive component.
#' @param sep_components Component separator.
#' @param ... Passed to `netmeta::netcomb()`.
#' @return An `nmaflow_fit` object.
#' @export
nma_component <- function(fit, inactive = NULL, sep_components = "+", ...) {
  if (!inherits(fit, "nmaflow_fit") || fit$engine != "netmeta")
    stop("The stable component-NMA implementation currently requires a netmeta-backed fit.", call. = FALSE)
  .nma_require("netmeta", "component NMA")
  cn <- netmeta::netcomb(fit$fit, inactive = inactive, sep.comps = sep_components, ...)
  .nma_new("nmaflow_fit", engine = "netmeta", framework = "frequentist_component",
           fit = cn, source_fit = fit, measure = fit$measure,
           specification = list(inactive = inactive, sep_components = sep_components))
}

#' Fit a Bayesian time-course Model-Based NMA
#'
#' @param data Arm-level data in MBNMAtime format.
#' @param fun A time-course function object accepted by `MBNMAtime::mb.run()`.
#' @param method Common/random setting recorded for the nmaFlow object.
#' @param reference Optional reference treatment for MBNMAtime when treatments are character labels.
#' @param ... Additional arguments to `MBNMAtime::mb.run()`.
#' @return An `nmaflow_fit` object.
#' @export
nma_time <- function(data, fun, method = "random", reference = NULL, ...) {
  .nma_require("MBNMAtime", "time-course NMA")
  method <- match.arg(method, c("common", "random"))
  if (is.null(reference) && "treatment" %in% names(data) && is.character(data[["treatment"]])) {
    reference <- unique(data[["treatment"]])[1L]
  }
  net <- MBNMAtime::mb.network(data, reference = reference)
  fit <- MBNMAtime::mb.run(net, fun = fun, ...)
  .nma_new("nmaflow_fit", engine = "MBNMAtime", framework = "bayesian_time",
           fit = fit, network = net, data = data, specification = list(method = method, reference = reference))
}

#' Create a multinma network from individual-level data
#'
#' @param data Individual participant/plot/plant-level data frame.
#' @param study,treatment Study and treatment column names.
#' @param outcome Outcome column.
#' @param outcome_type `continuous`, `binary`, or `poisson`.
#' @param exposure Optional exposure/time-at-risk column for Poisson outcomes.
#' @param reference Optional reference treatment.
#' @return A multinma `nma_data` object.
#' @export
nma_ipd <- function(data, study, treatment, outcome,
                    outcome_type = c("continuous", "binary", "poisson"),
                    exposure = NULL, reference = NULL) {
  .nma_require("multinma", "IPD NMA")
  .nma_assert_data(data); outcome_type <- match.arg(outcome_type)
  .nma_assert_cols(data, c(study, treatment, outcome, exposure))
  if (outcome_type == "continuous")
    return(multinma::set_ipd(data, study = study, trt = treatment, y = outcome, trt_ref = reference))
  if (outcome_type == "binary")
    return(multinma::set_ipd(data, study = study, trt = treatment, r = outcome, trt_ref = reference))
  if (is.null(exposure)) stop("Poisson IPD setup requires `exposure`.", call. = FALSE)
  multinma::set_ipd(data, study = study, trt = treatment, r = outcome, E = exposure, trt_ref = reference)
}

#' Combine IPD and aggregate evidence for ML-NMR preparation
#'
#' @param ipd_network A network created by [nma_ipd()] or `multinma::set_ipd()`.
#' @param agd_network Aggregate network created by [nma_make_multinma()] or a multinma set function.
#' @param reference Optional common reference treatment.
#' @param integration Optional function that accepts the combined network and returns a
#'   network with integration points added by `multinma::add_integration()`.
#' @return Combined multinma network, optionally with integration points.
#' @export
nma_mlnmr <- function(ipd_network, agd_network, reference = NULL, integration = NULL) {
  .nma_require("multinma", "ML-NMR")
  if (is.null(reference)) net <- multinma::combine_network(ipd_network, agd_network)
  else net <- multinma::combine_network(ipd_network, agd_network, trt_ref = reference)
  if (!is.null(integration)) {
    if (!is.function(integration)) stop("`integration` must be a function.", call. = FALSE)
    net <- integration(net)
  }
  net
}

#' Fit cross-design and cross-format NMA
#'
#' @param ... Arguments passed to `crossnma::crossnma.model()`.
#' @param run Whether to run JAGS after building the model.
#' @param run_args Named list passed to `crossnma::crossnma()` when `run=TRUE`.
#' @return Model specification or fitted `nmaflow_fit`.
#' @export
nma_crossdesign <- function(..., run = TRUE, run_args = list()) {
  .nma_require("crossnma", "cross-design/cross-format NMA")
  model_args <- list(...)
  required <- c("trt", "study", "outcome", "n", "design")
  if (any(!required %in% names(model_args))) {
    stop("`nma_crossdesign()` requires named arguments passed to crossnma::crossnma.model(), including trt, study, outcome, n and design.", call. = FALSE)
  }
  for (nm in intersect(c("trt", "study", "outcome", "n", "design", "se"), names(model_args))) {
    if (is.character(model_args[[nm]]) && length(model_args[[nm]]) == 1L) {
      model_args[[nm]] <- as.name(model_args[[nm]])
    }
  }
  model <- do.call(crossnma::crossnma.model, model_args)
  if (!isTRUE(run)) return(model)
  fit_args <- c(list(x = model), run_args)
  fit <- do.call(crossnma::crossnma, fit_args)
  .nma_new("nmaflow_fit", engine = "crossnma", framework = "bayesian_crossdesign",
           fit = fit, model = model, specification = list(call = match.call()))
}

#' Fit multivariate meta-analysis as an advanced NMA building block
#'
#' @param y Matrix/vector of correlated effect estimates.
#' @param S Within-study covariance specification.
#' @param formula Optional mixmeta formula. If omitted, `y ~ 1` is used.
#' @param data Optional data frame.
#' @param ... Passed to `mixmeta::mixmeta()`.
#' @return An `nmaflow_fit` object.
#' @export
nma_multivariate <- function(y, S, formula = NULL, data = NULL, ...) {
  .nma_require("mixmeta", "multivariate NMA building block")
  if (is.null(formula)) {
    if (is.null(data)) {
      data <- if (!is.null(dim(y))) data.frame(.nma_y = I(y)) else data.frame(.nma_y = y)
    }
    formula <- .nma_y ~ 1
  }
  fit <- mixmeta::mixmeta(formula, S = S, data = data, ...)
  .nma_new("nmaflow_fit", engine = "mixmeta", framework = "multivariate", fit = fit,
           specification = list(call = match.call()))
}

#' SEM/MASEM orchestration for network-derived effect structures
#'
#' @param engine `metaSEM` or `lavaan`.
#' @param ... Arguments passed to the selected engine function.
#' @param function_name Exported function name. Defaults to `meta()` for metaSEM and `sem()` for lavaan.
#' @return An `nmaflow_fit` object.
#' @export
nma_sem <- function(engine = c("metaSEM", "lavaan"), ..., function_name = NULL) {
  engine <- match.arg(engine)
  .nma_require(engine, "SEM-based evidence synthesis")
  function_name <- function_name %||% if (engine == "metaSEM") "meta" else "sem"
  if (!function_name %in% getNamespaceExports(engine))
    stop(sprintf("Function '%s' is not exported by installed package '%s'.", function_name, engine), call. = FALSE)
  fit <- do.call(getExportedValue(engine, function_name), list(...))
  .nma_new("nmaflow_fit", engine = engine, framework = "sem_nma", fit = fit,
           specification = list(function_name = function_name))
}
