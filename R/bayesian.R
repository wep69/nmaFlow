#' Create a multinma aggregate-data network from nmaFlow data
#'
#' @param x Arm-level `nmaflow_data`.
#' @param outcome_type `continuous` or `binomial`.
#' @param y,se Continuous outcome mean and SE columns.
#' @param event,total Binary event and denominator columns.
#' @param sample_size Optional arm sample-size column for continuous outcomes.
#' @param reference Optional network reference treatment.
#' @return A `multinma::nma_data` object.
#' @export
nma_make_multinma <- function(x, outcome_type = c("continuous", "binomial"),
                              y = NULL, se = NULL, event = NULL, total = NULL,
                              sample_size = NULL, reference = NULL) {
  .nma_require("multinma", "Bayesian NMA data setup")
  x <- .nma_unwrap(x); outcome_type <- match.arg(outcome_type)
  if (x$format != "arm") stop("multinma arm setup requires arm-level data.", call. = FALSE)
  d <- x$data; st <- x$mapping$study; tr <- x$mapping$treatment
  .nma_assert_cols(d, c(st, tr, y, se, event, total, sample_size))

  if (outcome_type == "continuous") {
    if (is.null(y) || is.null(se)) stop("Continuous multinma setup requires `y` and `se` columns.", call. = FALSE)
    if (is.null(sample_size)) {
      return(rlang::inject(multinma::set_agd_arm(
        data = d, study = !!rlang::sym(st), trt = !!rlang::sym(tr),
        y = !!rlang::sym(y), se = !!rlang::sym(se), trt_ref = reference
      )))
    }
    return(rlang::inject(multinma::set_agd_arm(
      data = d, study = !!rlang::sym(st), trt = !!rlang::sym(tr),
      y = !!rlang::sym(y), se = !!rlang::sym(se), sample_size = !!rlang::sym(sample_size),
      trt_ref = reference
    )))
  }

  if (is.null(event) || is.null(total)) stop("Binomial multinma setup requires `event` and `total` columns.", call. = FALSE)
  rlang::inject(multinma::set_agd_arm(
    data = d, study = !!rlang::sym(st), trt = !!rlang::sym(tr),
    r = !!rlang::sym(event), n = !!rlang::sym(total), trt_ref = reference
  ))
}

#' Construct a prior distribution for multinma
#'
#' @param distribution One of `normal`, `half_normal`, `student_t`, `half_student_t`, `cauchy`, `half_cauchy`, or `exponential`.
#' @param ... Parameters passed to the corresponding multinma prior constructor.
#' @return A multinma prior object.
#' @export
nma_prior <- function(distribution = c("normal", "half_normal", "student_t", "half_student_t",
                                       "cauchy", "half_cauchy", "exponential"), ...) {
  .nma_require("multinma", "Bayesian prior specification")
  distribution <- match.arg(distribution)
  fun <- getExportedValue("multinma", distribution)
  fun(...)
}

#' Fit Bayesian NMA or network meta-regression with multinma
#'
#' @param network A `multinma::nma_data` object.
#' @param consistency `consistency`, `ume`, or `nodesplit`.
#' @param trt_effects `fixed` or `random`.
#' @param regression Optional one-sided regression formula.
#' @param prior_trt,prior_het,prior_reg Optional multinma prior objects.
#' @param ... Additional arguments passed to `multinma::nma()` including Stan controls.
#' @return An `nmaflow_fit` object; node-splitting returns the corresponding wrapped result.
#' @export
nma_fit_bayes <- function(network, consistency = c("consistency", "ume", "nodesplit"),
                          trt_effects = c("random", "fixed"), regression = NULL,
                          prior_trt = NULL, prior_het = NULL, prior_reg = NULL, ...) {
  .nma_require("multinma", "Bayesian NMA")
  consistency <- match.arg(consistency); trt_effects <- match.arg(trt_effects)
  args <- list(network = network, consistency = consistency, trt_effects = trt_effects,
               regression = regression)
  if (!is.null(prior_trt)) args$prior_trt <- prior_trt
  if (!is.null(prior_het)) args$prior_het <- prior_het
  if (!is.null(prior_reg)) args$prior_reg <- prior_reg
  args <- c(args, list(...))
  fit <- do.call(multinma::nma, args)
  .nma_check_backend_result(fit, "multinma model")
  framework <- if (consistency == "nodesplit") "bayesian_nodesplit" else "bayesian"
  .nma_new("nmaflow_fit", engine = "multinma", framework = framework, fit = fit,
           network = network, specification = list(consistency = consistency, trt_effects = trt_effects,
                                                   regression = regression))
}

#' Fit Bayesian node-splitting models
#' @param network A multinma network.
#' @param ... Passed to [nma_fit_bayes()].
#' @export
nma_nodesplit_bayes <- function(network, ...) {
  nma_fit_bayes(network, consistency = "nodesplit", ...)
}

#' Posterior/model-fit checking for multinma models
#'
#' @param fit Bayesian `nmaflow_fit`.
#' @param type One of `dic`, `loo`, `waic`, or `summary`.
#' @param ... Backend arguments.
#' @export
nma_pp_check <- function(fit, type = c("dic", "loo", "waic", "summary"), ...) {
  if (!inherits(fit, "nmaflow_fit") || fit$engine != "multinma") stop("Requires a multinma-backed fit.", call. = FALSE)
  type <- match.arg(type)
  switch(type,
    dic = multinma::dic(fit$fit, ...),
    loo = {
      .nma_require("loo", "LOO-CV model comparison")
      loo::loo(fit$fit, ...)
    },
    waic = {
      .nma_require("loo", "WAIC model comparison")
      loo::waic(fit$fit, ...)
    },
    summary = summary(fit$fit, ...)
  )
}

#' Diagnose Bayesian sampling using posterior summaries
#'
#' @param fit Bayesian `nmaflow_fit`.
#' @param rhat_max Maximum acceptable R-hat used for the flag.
#' @param ess_min Minimum bulk/tail ESS used for the flag when available.
#' @return A list containing posterior summary and flagged rows.
#' @export
nma_mcmc_diagnose <- function(fit, rhat_max = 1.01, ess_min = 400) {
  if (!inherits(fit, "nmaflow_fit") || fit$engine != "multinma") stop("Requires a multinma-backed fit.", call. = FALSE)
  sm <- as.data.frame(summary(fit$fit))
  rhat_col <- intersect(c("Rhat", "rhat"), names(sm))
  bulk_col <- intersect(c("Bulk_ESS", "n_eff", "ESS_bulk"), names(sm))
  tail_col <- intersect(c("Tail_ESS", "ESS_tail"), names(sm))
  flag <- rep(FALSE, nrow(sm))
  if (length(rhat_col)) flag <- flag | (!is.na(sm[[rhat_col[1L]]]) & sm[[rhat_col[1L]]] > rhat_max)
  if (length(bulk_col)) flag <- flag | (!is.na(sm[[bulk_col[1L]]]) & sm[[bulk_col[1L]]] < ess_min)
  if (length(tail_col)) flag <- flag | (!is.na(sm[[tail_col[1L]]]) & sm[[tail_col[1L]]] < ess_min)
  list(summary = sm, flagged = sm[flag, , drop = FALSE], thresholds = c(rhat_max = rhat_max, ess_min = ess_min))
}

#' Bayesian relative-effect summary
#' @param fit Bayesian `nmaflow_fit`.
#' @param all_contrasts Return all treatment contrasts.
#' @param predictive Include the predictive distribution for random-effects models.
#' @param ... Additional arguments to `multinma::relative_effects()`.
#' @export
nma_bayes_summary <- function(fit, all_contrasts = TRUE, predictive = FALSE, ...) {
  if (!inherits(fit, "nmaflow_fit") || fit$engine != "multinma") stop("Requires a multinma-backed fit.", call. = FALSE)
  multinma::relative_effects(fit$fit, all_contrasts = all_contrasts,
                            predictive_distribution = predictive, ...)
}
