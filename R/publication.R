.nma_network_plot <- function(net) {
  .nma_require("ggplot2", "publication network plots")
  if (!inherits(net, "nmaflow_network")) stop("Expected nmaflow_network.", call. = FALSE)
  n <- nrow(net$nodes)
  if (!n) stop("Network has no nodes.", call. = FALSE)
  theta <- seq(0, 2 * pi, length.out = n + 1L)[-(n + 1L)]
  xy <- data.frame(treatment = net$nodes$treatment, x = cos(theta), y = sin(theta),
                   degree = net$nodes$weighted_degree, stringsAsFactors = FALSE)
  ed <- net$edges
  if (nrow(ed)) {
    ed$x <- xy$x[match(ed$trt1, xy$treatment)]; ed$y <- xy$y[match(ed$trt1, xy$treatment)]
    ed$xend <- xy$x[match(ed$trt2, xy$treatment)]; ed$yend <- xy$y[match(ed$trt2, xy$treatment)]
  }
  p <- ggplot2::ggplot()
  if (nrow(ed)) p <- p + ggplot2::geom_segment(data = ed, ggplot2::aes(x = x, y = y, xend = xend, yend = yend, linewidth = studies), alpha = 0.5)
  p <- p + ggplot2::geom_point(data = xy, ggplot2::aes(x = x, y = y, size = pmax(degree, 1))) +
    ggplot2::geom_text(data = xy, ggplot2::aes(x = x, y = y, label = treatment), nudge_y = 0.08, check_overlap = TRUE) +
    ggplot2::coord_equal() + ggplot2::theme_void() + ggplot2::labs(size = "Weighted degree", linewidth = "Studies")
  p
}

#' Create publication-oriented NMA figures
#'
#' @param x Data, network, fitted model, bootstrap, or ranking object.
#' @param type `network`, `forest`, `funnel`, `rank`, `inconsistency`, `bootstrap`, `dose`, or `time`.
#' @param file Optional output path. PDF/SVG are vector formats; TIFF/PNG use `dpi`.
#' @param width,height Figure dimensions in inches.
#' @param dpi Raster resolution.
#' @param ... Backend-specific plotting arguments.
#' @return Plot object invisibly when a backend draws directly.
#' @export
nma_plot <- function(x, type = c("network", "forest", "funnel", "rank", "inconsistency", "bootstrap", "dose", "time"),
                     file = NULL, width = 7, height = 5, dpi = 600, ...) {
  type <- match.arg(type)
  p <- NULL
  if (type == "network") {
    net <- if (inherits(x, "nmaflow_network")) x else if (inherits(x, "nmaflow_data")) nma_network(x) else NULL
    if (!is.null(net)) p <- .nma_network_plot(net)
    else if (inherits(x, "nmaflow_fit") && x$engine == "multinma") p <- plot(x$network, ...)
    else if (inherits(x, "nmaflow_fit") && x$engine == "netmeta") {
      p <- netmeta::netgraph(x$fit, ...); return(invisible(p))
    }
  } else if (type == "forest" && inherits(x, "nmaflow_fit")) {
    if (x$engine == "netmeta") {
      .nma_require("meta", "forest plots for netmeta models")
      return(invisible(meta::forest(x$fit, ...)))
    }
    if (x$engine == "NMA") {
      .nma_require("NMA", "forest plots for NMA models")
      return(invisible(NMA::nmaforest(x$fit, ...)))
    }
    if (x$engine == "multinma") {
      z <- multinma::relative_effects(x$fit, ...)
      return(plot(z))
    }
  } else if (type == "funnel" && inherits(x, "nmaflow_fit") && x$engine == "netmeta") {
    .nma_require("meta", "funnel plots for netmeta models")
    return(invisible(meta::funnel(x$fit, ...)))
  } else if (type == "rank" && inherits(x, "nmaflow_fit")) {
    p <- plot(nma_rank(x, ...))
  } else if (type == "inconsistency" && inherits(x, "nmaflow_fit")) {
    p <- nma_inconsistency(x, method = "heat", ...); return(invisible(p))
  } else if (type == "bootstrap" && inherits(x, "nmaflow_boot")) {
    .nma_require("ggplot2", "bootstrap plotting")
    d <- as.data.frame(x$replicates[x$success, , drop = FALSE])
    long <- stats::reshape(d, varying = names(d), v.names = "value", timevar = "term",
                           times = names(d), direction = "long")
    p <- ggplot2::ggplot(long, ggplot2::aes(x = value)) + ggplot2::geom_histogram(bins = 30) +
      ggplot2::facet_wrap(~ term, scales = "free") + ggplot2::theme_minimal(base_size = 11)
  } else if (type == "dose" && inherits(x, "nmaflow_fit")) {
    p <- plot(x$fit, ...)
  } else if (type == "time" && inherits(x, "nmaflow_fit")) {
    p <- plot(x$fit, ...)
  }
  if (is.null(p)) stop("Requested plot is unavailable for this object/backend.", call. = FALSE)
  if (!is.null(file)) {
    .nma_require("ggplot2", "figure export")
    ggplot2::ggsave(file, plot = p, width = width, height = height, dpi = dpi)
  }
  p
}

#' Create publication-oriented tables
#'
#' @param x nmaFlow object.
#' @param component `audit`, `issues`, `network`, `effects`, `bootstrap`, `rank`, `league`, or `summary`.
#' @param format `data.frame`, `gt`, or `flextable`.
#' @param digits Number of digits for numeric columns.
#' @param ... Passed to component extractors.
#' @return Table object.
#' @export
nma_table <- function(x, component = c("summary", "audit", "issues", "network", "effects", "bootstrap", "rank", "league"),
                      format = c("data.frame", "gt", "flextable"), digits = getOption("nmaFlow.digits", 3), ...) {
  component <- match.arg(component); format <- match.arg(format)
  tab <- switch(component,
    audit = nma_audit(x),
    issues = if (inherits(x, "nmaflow_validation")) x$issues else nma_validate(x, ...)$issues,
    network = { n <- if (inherits(x, "nmaflow_network")) x else nma_network(x); n$edges },
    effects = if (inherits(x, "nmaflow_effects")) x$data else stop("Use an nmaflow_effects object.", call. = FALSE),
    bootstrap = if (inherits(x, "nmaflow_boot")) x$ci else stop("Use an nmaflow_boot object.", call. = FALSE),
    rank = {
      r <- nma_rank(x, ...)
      if (inherits(r, "netrank")) {
        v <- r$Pscore.random %||% r$Pscore.common %||% r$Pscore.fixed
        if (is.character(v)) v <- r$Pscore.random
        if (is.null(v) || is.character(v)) stop("Could not extract numeric P-scores from the netmeta ranking object.", call. = FALSE)
        data.frame(treatment = names(v), P.score = as.numeric(v), stringsAsFactors = FALSE)
      } else {
        as.data.frame(r)
      }
    },
    league = as.data.frame(nma_league(x, ...)),
    summary = {
      if (inherits(x, "nmaflow_fit")) {
        s <- summary(x$fit); tryCatch(as.data.frame(s), error = function(e) data.frame(summary = utils::capture.output(print(s))))
      } else stop("Summary tables require a fitted model.", call. = FALSE)
    }
  )
  num <- vapply(tab, is.numeric, logical(1)); tab[num] <- lapply(tab[num], round, digits = digits)
  if (format == "data.frame") return(tab)
  if (format == "gt") { .nma_require("gt", "gt tables"); return(gt::gt(tab)) }
  .nma_require("flextable", "Word-ready tables"); flextable::flextable(tab)
}

#' Generate a league table
#' @param fit An `nmaflow_fit`.
#' @param ... Backend arguments.
#' @export
nma_league <- function(fit, ...) {
  if (!inherits(fit, "nmaflow_fit")) stop("Expected nmaflow_fit.", call. = FALSE)
  if (fit$engine %in% c("netmeta") && fit$framework == "frequentist_component") return(netmeta::netleague(fit$fit, ...))
  if (fit$engine == "netmeta") return(netmeta::netleague(fit$fit, ...))
  if (fit$engine == "NMA") return(NMA::nmaleague(fit$setup, ...))
  if (fit$engine == "multinma") return(multinma::relative_effects(fit$fit, all_contrasts = TRUE, ...))
  if (fit$engine == "crossnma") return(crossnma::league(fit$fit, ...))
  stop("League table unavailable for this backend.", call. = FALSE)
}

#' Summary-of-findings style evidence table
#' @param fit A netmeta-backed fit.
#' @param ... Passed to `netmeta::nettable()`.
#' @export
nma_sof <- function(fit, ...) {
  if (!inherits(fit, "nmaflow_fit") || fit$engine != "netmeta") stop("Current SoF table implementation requires netmeta.", call. = FALSE)
  netmeta::nettable(fit$fit, ...)
}

#' Write a reproducible Markdown analysis report skeleton
#'
#' @param x Data or fitted object.
#' @param file Output Markdown file.
#' @param title Report title.
#' @return Normalized output path invisibly.
#' @export
nma_report <- function(x, file = "nmaFlow_report.md", title = "Network Meta-Analysis") {
  lines <- c(
    paste0("# ", title), "",
    "## 1. Evidence question and estimand", "",
    "Document population/context, interventions, comparator policy, outcome, timepoint and target estimand.", "",
    "## 2. Data architecture and node policy", "",
    "Report arm/contrast format, coarse and fine nodes, units, timepoints, populations, duplicate handling and all reviewed mappings.", "",
    "## 3. Network preflight", "",
    "Report connected components, multi-arm studies, weakly supported edges, unavailable dispersion and go/no-go decision.", "",
    "## 4. Transitivity", "",
    "Describe prespecified effect modifiers and their distributions across comparisons.", "",
    "## 5. Statistical model", "",
    "Report framework, backend/version, effect measure, common/random effects, heterogeneity method or prior, and multi-arm handling.", "",
    "## 6. Inconsistency and diagnostics", "",
    "Report global/local inconsistency, node splitting where identifiable, convergence or numerical diagnostics, and predictive/model-fit checks.", "",
    "## 7. Sensitivity and bootstrap", "",
    "Report all alternative node, dispersion, zero-event, missing-data, outlier and bootstrap scenarios with seeds and replicate counts.", "",
    "## 8. Results", "",
    "Present effect estimates with uncertainty and prediction intervals where appropriate. Treat rankings as secondary to effect magnitude and uncertainty.", "",
    "## 9. Publication outputs", "",
    "Include network graph, forest/league table, inconsistency diagnostics, ranking uncertainty and dose/time curves where applicable.", "",
    "## 10. Reproducibility", "",
    paste0("Generated by nmaFlow at ", format(Sys.time(), tz = "UTC", usetz = TRUE), "."), "",
    "Append `nma_session()` output and the data audit log."
  )
  writeLines(lines, file)
  invisible(normalizePath(file, mustWork = FALSE))
}
