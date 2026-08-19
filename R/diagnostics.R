#' Assess transitivity using effect-modifier distributions
#'
#' @param x An `nmaflow_data` object.
#' @param modifiers Character vector of effect-modifier columns.
#' @param by Treatment node representation.
#' @return An object with treatment-specific summaries and overlap flags.
#' @export
nma_transitivity <- function(x, modifiers, by = c("treatment", "coarse", "fine")) {
  x <- .nma_unwrap(x); by <- match.arg(by)
  col <- switch(by, treatment = x$mapping$treatment, coarse = x$mapping$node_coarse, fine = x$mapping$node_fine)
  if (is.null(col)) stop("Requested node representation is unavailable.", call. = FALSE)
  .nma_assert_cols(x$data, c(col, modifiers))
  out <- list()
  for (v in modifiers) {
    z <- x$data[[v]]; tr <- as.character(x$data[[col]])
    if (is.numeric(z)) {
      sp <- split(z, tr)
      out[[v]] <- data.frame(
        treatment = names(sp), n = vapply(sp, function(a) sum(is.finite(a)), integer(1)),
        mean = vapply(sp, function(a) mean(a, na.rm = TRUE), numeric(1)),
        sd = vapply(sp, function(a) stats::sd(a, na.rm = TRUE), numeric(1)),
        min = vapply(sp, function(a) min(a, na.rm = TRUE), numeric(1)),
        max = vapply(sp, function(a) max(a, na.rm = TRUE), numeric(1)), stringsAsFactors = FALSE
      )
    } else {
      tab <- table(tr, as.character(z), useNA = "ifany")
      out[[v]] <- prop.table(tab, 1L)
    }
  }
  .nma_new("nmaflow_transitivity", modifiers = modifiers, summaries = out,
           note = "Transitivity is a scientific assumption. Distributional overlap is diagnostic evidence, not a formal proof.")
}

#' Assess network inconsistency with the fitted backend
#'
#' @param fit An `nmaflow_fit`.
#' @param method `netsplit`, `design`, or `heat` for netmeta; `global` or `local` for NMA.
#' @param ... Backend arguments.
#' @export
nma_inconsistency <- function(fit, method = NULL, ...) {
  if (!inherits(fit, "nmaflow_fit")) stop("Expected `nmaflow_fit`.", call. = FALSE)
  if (fit$engine == "netmeta") {
    method <- method %||% "netsplit"
    return(switch(method,
      netsplit = netmeta::netsplit(fit$fit, ...),
      design = netmeta::decomp.design(fit$fit, ...),
      heat = netmeta::netheat(fit$fit, ...),
      stop("Unknown netmeta inconsistency method.", call. = FALSE)
    ))
  }
  if (fit$engine == "NMA") {
    method <- method %||% "global"
    return(switch(method,
      global = NMA::global.ict(fit$setup, ...),
      local = NMA::local.ict(fit$setup, ...),
      stop("Unknown NMA inconsistency method.", call. = FALSE)
    ))
  }
  if (fit$engine == "multinma") return(nma_nodesplit_bayes(fit$network, ...))
  stop("Inconsistency method unavailable for this backend.", call. = FALSE)
}

#' Node-splitting convenience wrapper
#' @param fit An `nmaflow_fit`.
#' @param ... Backend arguments.
#' @export
nma_nodesplit <- function(fit, ...) {
  if (fit$engine == "netmeta") return(netmeta::netsplit(fit$fit, ...))
  if (fit$engine == "multinma") return(nma_nodesplit_bayes(fit$network, ...))
  if (fit$engine == "NMA") return(NMA::sidesplit(fit$setup, ...))
  stop("Node splitting unavailable for this backend.", call. = FALSE)
}

#' Design-by-treatment inconsistency decomposition
#' @param fit A netmeta-backed fit.
#' @param ... Backend arguments.
#' @export
nma_design_by_treatment <- function(fit, ...) {
  if (!inherits(fit, "nmaflow_fit") || fit$engine != "netmeta") stop("Requires a netmeta-backed fit.", call. = FALSE)
  netmeta::decomp.design(fit$fit, ...)
}

#' Net heat plot
#' @param fit A netmeta-backed fit.
#' @param ... Backend arguments.
#' @export
nma_netheat <- function(fit, ...) {
  if (!inherits(fit, "nmaflow_fit") || fit$engine != "netmeta") stop("Requires a netmeta-backed fit.", call. = FALSE)
  netmeta::netheat(fit$fit, ...)
}

#' Identify closed loops in the observed network
#'
#' @param x `nmaflow_data` or `nmaflow_network`.
#' @param max_length Maximum loop length searched; default 4.
#' @return Data frame of unique simple cycles found by a bounded depth-first search.
#' @export
nma_loop_check <- function(x, max_length = 4L) {
  net <- if (inherits(x, "nmaflow_network")) x else nma_network(x)
  nodes <- net$nodes$treatment; es <- net$edges
  adj <- stats::setNames(vector("list", length(nodes)), nodes)
  if (nrow(es)) for (i in seq_len(nrow(es))) {
    adj[[es$trt1[i]]] <- unique(c(adj[[es$trt1[i]]], es$trt2[i]))
    adj[[es$trt2[i]]] <- unique(c(adj[[es$trt2[i]]], es$trt1[i]))
  }
  cycles <- character()
  canon <- function(path) {
    z <- path[-length(path)]
    rots <- lapply(seq_along(z), function(i) c(z[i:length(z)], z[seq_len(i - 1L)]))
    revz <- rev(z); rots2 <- lapply(seq_along(revz), function(i) c(revz[i:length(revz)], revz[seq_len(i - 1L)]))
    min(vapply(c(rots, rots2), paste, collapse = " -> ", FUN.VALUE = character(1)))
  }
  dfs <- function(start, current, path) {
    if (length(path) > max_length) return()
    for (nxt in adj[[current]]) {
      if (nxt == start && length(path) >= 3L) cycles <<- unique(c(cycles, canon(c(path, start))))
      else if (!nxt %in% path && length(path) < max_length) dfs(start, nxt, c(path, nxt))
    }
  }
  for (s in nodes) dfs(s, s, s)
  data.frame(loop = sort(cycles), stringsAsFactors = FALSE)
}
