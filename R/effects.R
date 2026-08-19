#' Recommend an effect measure from outcome-scale properties
#'
#' @param scale_type One of `ratio`, `interval`, `ordinal`, `count`, `proportion`, `binary`, or `time_to_event`.
#' @param can_be_negative Logical indicating whether observed values can be negative.
#' @param same_scale Logical; for continuous outcomes, whether all studies use a common interpretable scale.
#' @return A character recommendation.
#' @export
nma_measure <- function(scale_type, can_be_negative = FALSE, same_scale = TRUE) {
  scale_type <- match.arg(scale_type, c("ratio", "interval", "ordinal", "count", "proportion", "binary", "time_to_event"))
  switch(scale_type,
    ratio = if (isTRUE(can_be_negative)) "MD" else "ROM",
    interval = if (isTRUE(same_scale)) "MD" else "SMD",
    ordinal = "OR",
    count = "IRR",
    proportion = "RR_or_OR",
    binary = "RR_or_OR",
    time_to_event = "HR"
  )
}

#' Harmonize effect direction
#'
#' @param x Numeric outcome or effect vector.
#' @param higher_is_better Logical indicating the original favorable direction.
#' @param target Either `higher` or `lower`.
#' @return Direction-harmonized numeric vector.
#' @export
nma_direction <- function(x, higher_is_better = TRUE, target = c("higher", "lower")) {
  target <- match.arg(target)
  flip <- xor(isTRUE(higher_is_better), identical(target, "higher"))
  if (flip) -x else x
}

#' Recover or document dispersion from reported summary information
#'
#' @param n Sample size.
#' @param mean Optional mean.
#' @param sd Reported SD when available.
#' @param se Reported SE when available.
#' @param ci_lower,ci_upper Confidence limits for the mean.
#' @param conf_level Confidence level associated with the interval.
#' @param mse Optional ANOVA mean square error.
#' @param cv_percent Optional experimental coefficient of variation reported by the study.
#' @param lsd Optional Fisher LSD value.
#' @param hsd Optional Tukey HSD value.
#' @param df_error Error degrees of freedom required for LSD recovery.
#' @param studentized_q Tukey studentized range critical value required for HSD recovery.
#' @return A list with effective SD, method, approximation flag and formula note.
#' @export
nma_dispersion <- function(n, mean = NA_real_, sd = NA_real_, se = NA_real_,
                           ci_lower = NA_real_, ci_upper = NA_real_, conf_level = 0.95,
                           mse = NA_real_, cv_percent = NA_real_, lsd = NA_real_, hsd = NA_real_,
                           df_error = NA_real_, studentized_q = NA_real_) {
  if (!is.finite(n) || n <= 0) stop("`n` must be positive.", call. = FALSE)
  if (is.finite(sd) && sd > 0) return(list(sd = sd, method = "reported_sd", approximate = FALSE, formula = "reported"))
  if (is.finite(se) && se > 0) return(list(sd = se * sqrt(n), method = "se_to_sd", approximate = FALSE, formula = "SD = SE * sqrt(n)"))
  if (is.finite(ci_lower) && is.finite(ci_upper) && ci_upper > ci_lower) {
    z <- stats::qnorm(1 - (1 - conf_level) / 2)
    se2 <- (ci_upper - ci_lower) / (2 * z)
    return(list(sd = se2 * sqrt(n), method = "ci_to_sd", approximate = TRUE,
                formula = "SE = (upper-lower)/(2*z); SD = SE*sqrt(n)"))
  }
  if (is.finite(mse) && mse > 0) return(list(sd = sqrt(mse), method = "mse_to_sd", approximate = TRUE, formula = "SD = sqrt(MSE)"))
  if (is.finite(cv_percent) && cv_percent > 0 && is.finite(mean)) {
    return(list(sd = abs(mean) * cv_percent / 100, method = "cv_to_sd", approximate = TRUE,
                formula = "SD = |mean| * CV/100; valid only for article-reported experimental CV"))
  }
  if (is.finite(lsd) && lsd > 0 && is.finite(df_error) && df_error > 0) {
    tcrit <- stats::qt(0.975, df = df_error)
    mse2 <- (lsd / tcrit)^2 * n / 2
    return(list(sd = sqrt(mse2), method = "lsd_to_sd", approximate = TRUE,
                formula = "MSE=(LSD/t)^2*n/2; SD=sqrt(MSE)"))
  }
  if (is.finite(hsd) && hsd > 0 && is.finite(studentized_q) && studentized_q > 0) {
    mse2 <- (hsd / studentized_q)^2 * n
    return(list(sd = sqrt(mse2), method = "hsd_to_sd", approximate = TRUE,
                formula = "MSE=(HSD/q)^2*n; SD=sqrt(MSE)"))
  }
  list(sd = NA_real_, method = "unavailable", approximate = NA, formula = "No defensible recovery source supplied")
}

#' Check multi-arm evidence structure
#'
#' @param x An `nmaflow_data` object.
#' @return Data frame with arm counts and flags by study.
#' @export
nma_multiarm_check <- function(x) {
  x <- .nma_unwrap(x)
  if (x$format != "arm") stop("Multi-arm checks require arm-level data.", call. = FALSE)
  d <- x$data; st <- x$mapping$study; tr <- x$mapping$treatment
  sp <- split(as.character(d[[tr]]), as.character(d[[st]]))
  data.frame(
    study = names(sp),
    n_rows = vapply(sp, length, integer(1)),
    n_treatments = vapply(sp, function(z) length(unique(z[!is.na(z)])), integer(1)),
    multi_arm = vapply(sp, function(z) length(unique(z[!is.na(z)])) > 2L, logical(1)),
    repeated_node = vapply(sp, function(z) any(duplicated(z[!is.na(z)])), logical(1)),
    stringsAsFactors = FALSE
  )
}

#' Convert arm-level summaries to contrast-level effect estimates
#'
#' This helper is intentionally transparent and is best suited to teaching, diagnostics,
#' bootstrap workflows, and backends that accept contrast-level inputs. Multi-arm studies
#' are preserved by `study` and are flagged because pairwise contrasts are correlated.
#'
#' @param x An arm-level `nmaflow_data` object.
#' @param measure One of `MD`, `ROM`, `RR`, `OR`, or `RD`.
#' @param mean,sd,n Columns for continuous outcomes.
#' @param event,total Columns for binary outcomes.
#' @param continuity Continuity correction used for zero-cell log RR/OR calculations.
#' @return An `nmaflow_effects` object.
#' @export
nma_effects <- function(x, measure = c("MD", "ROM", "RR", "OR", "RD"),
                        mean = NULL, sd = NULL, n = NULL, event = NULL, total = NULL,
                        continuity = 0.5) {
  x <- .nma_unwrap(x); measure <- match.arg(measure)
  if (x$format != "arm") stop("nma_effects() currently requires arm-level data.", call. = FALSE)
  d <- x$data; m <- x$mapping
  .nma_assert_cols(d, c(mean, sd, n, event, total))
  continuous <- measure %in% c("MD", "ROM")
  if (continuous && any(vapply(list(mean, sd, n), is.null, logical(1)))) stop("Continuous effects require `mean`, `sd`, and `n`.", call. = FALSE)
  if (!continuous && any(vapply(list(event, total), is.null, logical(1)))) stop("Binary effects require `event` and `total`.", call. = FALSE)

  sp <- split(seq_len(nrow(d)), as.character(d[[m$study]]))
  out <- list(); k <- 0L
  for (s in names(sp)) {
    idx <- sp[[s]]; z <- d[idx, , drop = FALSE]
    if (nrow(z) < 2L) next
    cmb <- utils::combn(seq_len(nrow(z)), 2L)
    for (j in seq_len(ncol(cmb))) {
      i1 <- cmb[1L, j]; i2 <- cmb[2L, j]
      t1 <- as.character(z[[m$treatment]][i1]); t2 <- as.character(z[[m$treatment]][i2])
      if (identical(t1, t2)) next
      k <- k + 1L
      if (measure == "MD") {
        te <- z[[mean]][i2] - z[[mean]][i1]
        se <- sqrt(z[[sd]][i1]^2 / z[[n]][i1] + z[[sd]][i2]^2 / z[[n]][i2])
      } else if (measure == "ROM") {
        if (z[[mean]][i1] <= 0 || z[[mean]][i2] <= 0) { te <- se <- NA_real_ } else {
          te <- log(z[[mean]][i2] / z[[mean]][i1])
          se <- sqrt(z[[sd]][i1]^2 / (z[[n]][i1] * z[[mean]][i1]^2) +
                     z[[sd]][i2]^2 / (z[[n]][i2] * z[[mean]][i2]^2))
        }
      } else {
        a <- z[[event]][i2]; n2 <- z[[total]][i2]; c <- z[[event]][i1]; n1 <- z[[total]][i1]
        b <- n2 - a; dd <- n1 - c
        if (measure %in% c("RR", "OR") && any(c(a, b, c, dd) == 0)) {
          a <- a + continuity; b <- b + continuity; c <- c + continuity; dd <- dd + continuity
          n2 <- a + b; n1 <- c + dd
        }
        if (measure == "RR") {
          te <- log((a / n2) / (c / n1))
          se <- sqrt(1 / a - 1 / n2 + 1 / c - 1 / n1)
        } else if (measure == "OR") {
          te <- log((a / b) / (c / dd))
          se <- sqrt(1 / a + 1 / b + 1 / c + 1 / dd)
        } else {
          te <- a / n2 - c / n1
          se <- sqrt((a / n2) * (1 - a / n2) / n2 + (c / n1) * (1 - c / n1) / n1)
        }
      }
      out[[k]] <- data.frame(study = s, treat1 = t1, treat2 = t2, TE = te, seTE = se,
                             measure = measure, multi_arm = nrow(z) > 2L, stringsAsFactors = FALSE)
    }
  }
  dat <- if (length(out)) do.call(rbind, out) else data.frame()
  .nma_new("nmaflow_effects", data = dat, measure = measure,
           note = "Contrasts from the same multi-arm study are correlated; use a backend that models this dependence.")
}

#' @export
print.nmaflow_effects <- function(x, ...) {
  cat("nmaFlow contrast data\n")
  cat("  measure:", x$measure, "\n")
  cat("  contrasts:", nrow(x$data), "\n")
  if (nrow(x$data)) print(utils::head(x$data, 10L), row.names = FALSE)
  invisible(x)
}

#' Summarize effect-measure readiness
#'
#' @param x An `nmaflow_data` object.
#' @param ... Passed to [nma_validate()].
#' @return A list containing validation, geometry, multi-arm summary and go/no-go status.
#' @export
nma_readiness <- function(x, ...) {
  pre <- nma_preflight(x, ...)
  ma <- if (x$format == "arm") nma_multiarm_check(x) else NULL
  .nma_new("nmaflow_readiness", status = pre$status, validation = pre$validation,
           network = pre$network, multiarm = ma)
}
