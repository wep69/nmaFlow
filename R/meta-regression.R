#' Fit network meta-regression through the active backend
#'
#' @param fit An `nmaflow_fit`.
#' @param regression For multinma, a one-sided formula. For other backends pass the
#'   backend-specific covariate specification in `...`.
#' @param ... Additional backend arguments.
#' @return Backend model, wrapped when refitting multinma.
#' @export
nma_regression <- function(fit, regression = NULL, ...) {
  if (!inherits(fit, "nmaflow_fit")) stop("Expected `nmaflow_fit`.", call. = FALSE)
  if (fit$engine == "netmeta") return(netmeta::netmetareg(fit$fit, ...))
  if (fit$engine == "NMA") return(NMA::nmareg(fit$setup, ...))
  if (fit$engine == "multinma") {
    return(nma_fit_bayes(fit$network,
                         consistency = fit$specification$consistency %||% "consistency",
                         trt_effects = fit$specification$trt_effects %||% "random",
                         regression = regression, ...))
  }
  stop("Network meta-regression is not implemented for this backend.", call. = FALSE)
}

#' Summarize candidate effect modifiers before meta-regression
#'
#' @param x An `nmaflow_data` object.
#' @param modifiers Candidate effect modifiers.
#' @return The transitivity summary plus missingness by modifier.
#' @export
nma_effect_modifier <- function(x, modifiers) {
  x <- .nma_unwrap(x); .nma_assert_cols(x$data, modifiers)
  miss <- data.frame(
    modifier = modifiers,
    missing_n = vapply(modifiers, function(v) sum(is.na(x$data[[v]])), integer(1)),
    missing_fraction = vapply(modifiers, function(v) mean(is.na(x$data[[v]])), numeric(1)),
    stringsAsFactors = FALSE
  )
  list(transitivity = nma_transitivity(x, modifiers), missingness = miss)
}

#' Plot a network meta-regression relationship
#'
#' @param data Data frame containing observed study-level covariate and effect estimates.
#' @param x Covariate column.
#' @param y Effect-estimate column.
#' @param se Optional SE column used to display uncertainty.
#' @param treatment Optional comparison/treatment column for grouping.
#' @return A ggplot object.
#' @export
nma_regression_plot <- function(data, x, y, se = NULL, treatment = NULL) {
  .nma_require("ggplot2", "network meta-regression plotting")
  .nma_assert_data(data); .nma_assert_cols(data, c(x, y, se, treatment))
  p <- ggplot2::ggplot(data, ggplot2::aes_string(x = x, y = y)) +
    ggplot2::geom_point() + ggplot2::theme_minimal(base_size = 11) +
    ggplot2::labs(x = x, y = y)
  if (!is.null(se)) {
    data$.nma_ymin <- data[[y]] - 1.96 * data[[se]]
    data$.nma_ymax <- data[[y]] + 1.96 * data[[se]]
    p <- ggplot2::ggplot(data, ggplot2::aes_string(x = x, y = y)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = .nma_ymin, ymax = .nma_ymax), width = 0) +
      ggplot2::geom_point() + ggplot2::theme_minimal(base_size = 11) + ggplot2::labs(x = x, y = y)
  }
  if (!is.null(treatment)) p <- p + ggplot2::aes_string(shape = treatment)
  p
}

#' Compare a set of network meta-regression models
#' @param fits Named list of backend models or nmaFlow fits.
#' @return Data frame of AIC/BIC when available.
#' @export
nma_regression_compare <- function(fits) {
  if (!is.list(fits) || !length(fits)) stop("Provide a non-empty list of models.", call. = FALSE)
  if (is.null(names(fits))) names(fits) <- paste0("model", seq_along(fits))
  do.call(rbind, lapply(names(fits), function(nm) {
    z <- fits[[nm]]; if (inherits(z, "nmaflow_fit")) z <- z$fit
    data.frame(model = nm,
               AIC = tryCatch(stats::AIC(z), error = function(e) NA_real_),
               BIC = tryCatch(stats::BIC(z), error = function(e) NA_real_),
               stringsAsFactors = FALSE)
  }))
}
