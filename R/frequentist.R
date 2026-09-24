.nma_effect_data <- function(x) {
  if (inherits(x, "nmaflow_effects")) return(x$data)
  if (is.data.frame(x)) return(x)
  stop("Provide an `nmaflow_effects` object or contrast data frame.", call. = FALSE)
}

#' Fit frequentist network meta-analysis with netmeta
#'
#' @param x Contrast data from [nma_effects()] or a data frame.
#' @param studlab,treat1,treat2,TE,seTE Column names.
#' @param sm Summary measure passed to `netmeta`.
#' @param common,random Fit common- and/or random-effects models.
#' @param reference Reference treatment.
#' @param method_tau Heterogeneity estimator.
#' @param ... Additional arguments passed to `netmeta::netmeta()`.
#' @return An `nmaflow_fit` object.
#' @export
nma_fit_freq <- function(x, studlab = "study", treat1 = "treat1", treat2 = "treat2",
                         TE = "TE", seTE = "seTE", sm = NULL, common = FALSE, random = TRUE,
                         reference = NULL, method_tau = "REML", ...) {
  .nma_require("netmeta", "frequentist NMA")
  d <- .nma_effect_data(x)
  .nma_assert_cols(d, c(studlab, treat1, treat2, TE, seTE))
  sm <- sm %||% if (inherits(x, "nmaflow_effects")) x$measure else "MD"
  netmeta_args <- list(
    TE = d[[TE]], seTE = d[[seTE]], treat1 = d[[treat1]], treat2 = d[[treat2]],
    studlab = d[[studlab]], data = d, sm = sm, common = common, random = random,
    method.tau = method_tau, ...
  )
  if (!is.null(reference)) netmeta_args$reference.group <- reference
  fit <- do.call(netmeta::netmeta, netmeta_args)
  .nma_new("nmaflow_fit", engine = "netmeta", framework = "frequentist", fit = fit,
           data = d, measure = sm, specification = list(common = common, random = random,
                                                        reference = reference, method_tau = method_tau))
}

#' Fit NMA using the NMA package improved REML framework
#'
#' @param x Arm-level `nmaflow_data` object.
#' @param measure `OR`, `RR`, `RD`, `MD`, `SMD`, `HR`, or `SPD`.
#' @param ref Reference treatment.
#' @param event,total,mean,sd,n Column names as required by the measure.
#' @param covariates Optional covariate column names for later meta-regression/transitivity.
#' @param method `NH`, `REML`, or `fixed`.
#' @param eform Exponentiate results where appropriate.
#' @return An `nmaflow_fit` object containing both setup and fit objects.
#' @export
nma_fit_nma <- function(x, measure = c("MD", "SMD", "OR", "RR", "RD", "HR", "SPD"), ref,
                        event = NULL, total = NULL, mean = NULL, sd = NULL, n = NULL,
                        covariates = NULL, method = c("NH", "REML", "fixed"), eform = FALSE) {
  .nma_require("NMA", "Noma-Hamura multivariate NMA")
  x <- .nma_unwrap(x); measure <- match.arg(measure); method <- match.arg(method)
  if (x$format != "arm") stop("NMA::setup requires arm-level data.", call. = FALSE)
  d <- x$data; st <- x$mapping$study; tr <- x$mapping$treatment
  .nma_assert_cols(d, c(st, tr, event, total, mean, sd, n, covariates))
  args <- list(
    study = as.name(st), trt = as.name(tr), measure = measure,
    ref = ref, data = d
  )
  if (!is.null(event)) args$d <- as.name(event)
  if (!is.null(total)) args$n <- as.name(total)
  if (!is.null(mean)) args$m <- as.name(mean)
  if (!is.null(sd)) args$s <- as.name(sd)
  if (!is.null(n)) args$n <- as.name(n)
  if (!is.null(covariates)) {
    args$z <- as.call(c(list(as.name("c")), lapply(covariates, function(v) v)))
  }
  setup <- do.call(NMA::setup, args)
  fit <- NMA::nma(setup, eform = eform, method = method)
  .nma_new("nmaflow_fit", engine = "NMA", framework = "frequentist", fit = fit,
           setup = setup, data = d, measure = measure, specification = list(method = method, ref = ref))
}

#' Fit additive component NMA using netmeta
#'
#' @param fit A standard netmeta-backed `nmaflow_fit`.
#' @param inactive Optional inactive component.
#' @param sep_components Component separator.
#' @param ... Passed to `netmeta::netcomb()`.
#' @return An `nmaflow_fit` component model.
#' @export
nma_fit_additive <- function(fit, inactive = NULL, sep_components = "+", ...) {
  nma_component(fit, inactive = inactive, sep_components = sep_components, ...)
}

#' Extract direct evidence when supported by the backend
#' @param fit An `nmaflow_fit` object.
#' @param ... Backend arguments.
#' @export
nma_direct <- function(fit, ...) {
  if (!inherits(fit, "nmaflow_fit")) stop("Expected `nmaflow_fit`.", call. = FALSE)
  if (fit$engine == "netmeta") return(netmeta::netpairwise(fit$fit, ...))
  stop("Direct-evidence extractor is not implemented for this backend.", call. = FALSE)
}

#' Contribution matrix or weights
#' @param fit An `nmaflow_fit` object.
#' @param ... Backend arguments.
#' @export
nma_contribution <- function(fit, ...) {
  if (!inherits(fit, "nmaflow_fit")) stop("Expected `nmaflow_fit`.", call. = FALSE)
  if (fit$engine == "netmeta") {
    if (exists("netcontrib", envir = asNamespace("netmeta"), inherits = FALSE))
      return(get("netcontrib", envir = asNamespace("netmeta"))(fit$fit, ...))
    return(fit$fit$H %||% NULL)
  }
  if (fit$engine == "NMA") return(NMA::nmaweight(fit$setup, ...))
  stop("Contribution analysis is not implemented for this backend.", call. = FALSE)
}

#' Predict treatment effects or intervals using the fitted backend
#' @param fit An `nmaflow_fit` object.
#' @param ... Passed to the backend prediction method.
#' @export
nma_predict <- function(fit, ...) {
  if (!inherits(fit, "nmaflow_fit")) stop("Expected `nmaflow_fit`.", call. = FALSE)
  if (fit$engine == "netmeta") {
    stop("Prediction is not available for netmeta fits. Use `nma_league()`, `nma_inconsistency()` or a backend with prediction support.", call. = FALSE)
  }
  stats::predict(fit$fit, ...)
}

#' Frequentist model summary
#' @param fit An `nmaflow_fit` object.
#' @param ... Passed to `summary()`.
#' @export
nma_freq_summary <- function(fit, ...) {
  if (!inherits(fit, "nmaflow_fit") || fit$framework != "frequentist")
    stop("Expected a frequentist `nmaflow_fit`.", call. = FALSE)
  summary(fit$fit, ...)
}

#' @export
print.nmaflow_fit <- function(x, ...) {
  cat("nmaFlow fit\n")
  cat("  engine:", x$engine, "\n")
  cat("  framework:", x$framework, "\n")
  if (!is.null(x$measure)) cat("  measure:", x$measure, "\n")
  print(x$fit)
  invisible(x)
}

#' @export
summary.nmaflow_fit <- function(object, ...) summary(object$fit, ...)
