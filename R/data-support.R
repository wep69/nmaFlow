#' Apply an explicit unit conversion with an audit trail
#'
#' `nma_convert_units()` never guesses unit conversions. The analyst supplies the
#' source unit, target unit, and conversion rule. This is deliberately conservative
#' because an incorrect automatic conversion can create a statistically coherent but
#' scientifically invalid network.
#'
#' @param x An `nmaflow_data` object.
#' @param value Column containing the numeric value to convert.
#' @param unit Column containing units.
#' @param from,to Source and target unit labels.
#' @param multiplier Multiplicative conversion factor.
#' @param offset Additive offset applied after multiplication.
#' @param reviewed Must be `TRUE`; conversion is an explicit scientific decision.
#' @param reason Text documenting the conversion source or rationale.
#' @return Updated `nmaflow_data` object with audit trail.
#' @export
nma_convert_units <- function(x, value, unit, from, to, multiplier, offset = 0,
                              reviewed = FALSE, reason = NULL) {
  x <- .nma_unwrap(x)
  if (!isTRUE(reviewed)) {
    stop("Unit conversion requires reviewed = TRUE and an explicit rule.", call. = FALSE)
  }
  .nma_assert_cols(x$data, c(value, unit))
  if (!is.numeric(multiplier) || length(multiplier) != 1L || !is.finite(multiplier)) {
    stop("multiplier must be one finite numeric value.", call. = FALSE)
  }
  if (!is.numeric(offset) || length(offset) != 1L || !is.finite(offset)) {
    stop("offset must be one finite numeric value.", call. = FALSE)
  }
  idx <- !is.na(x$data[[unit]]) & as.character(x$data[[unit]]) == from
  if (!any(idx)) {
    warning("No rows matched unit = '", from, "'.", call. = FALSE)
    return(x)
  }
  old <- x$data[[value]][idx]
  if (!is.numeric(old)) stop("The value column must be numeric for unit conversion.", call. = FALSE)
  x$data[[value]][idx] <- old * multiplier + offset
  x$data[[unit]][idx] <- to
  x$audit <- rbind(
    x$audit,
    .nma_audit_row(
      step = "repair",
      action = "unit_conversion",
      field = paste(value, unit, sep = "/"),
      before = paste0(from, "; n=", sum(idx)),
      after = paste0(to, "; multiplier=", multiplier, "; offset=", offset),
      rows = sum(idx),
      reason = reason %||% "Explicit reviewed unit conversion",
      automatic = FALSE
    )
  )
  x
}

#' Apply a reviewed fine/coarse treatment-node policy
#'
#' A node dictionary should contain one row per recognized treatment label and may
#' include both a fine-grained node and a coarse synthesis node. No fuzzy matching is
#' applied here: labels must match the dictionary exactly after optional conservative
#' text normalization performed by `nma_repair()`.
#'
#' @param x An `nmaflow_data` object.
#' @param dictionary Data frame containing node definitions.
#' @param label_col First/arm-level treatment label column in `x$data`.
#' @param label_col2 Optional second treatment label for contrast-level evidence; defaults to `treatment2`.
#' @param dictionary_label Label column in `dictionary`.
#' @param fine_col,coarse_col Columns in `dictionary` containing fine and coarse nodes.
#' @param reviewed Must be `TRUE`.
#' @param reason Audit rationale.
#' @return Updated `nmaflow_data` object.
#' @export
nma_node_policy <- function(x, dictionary, label_col = NULL, label_col2 = NULL,
                            dictionary_label = "label",
                            fine_col = "fine", coarse_col = "coarse",
                            reviewed = FALSE, reason = NULL) {
  x <- .nma_unwrap(x)
  if (!isTRUE(reviewed)) stop("Node policy requires reviewed = TRUE.", call. = FALSE)
  if (!is.data.frame(dictionary)) stop("dictionary must be a data frame.", call. = FALSE)
  label_col <- label_col %||% x$mapping$treatment
  .nma_assert_cols(x$data, label_col)
  .nma_assert_cols(dictionary, c(dictionary_label, fine_col, coarse_col))

  key <- match(as.character(x$data[[label_col]]), as.character(dictionary[[dictionary_label]]))
  unmatched <- is.na(key) & !is.na(x$data[[label_col]])
  if (any(unmatched)) {
    warning(sum(unmatched), " treatment labels were not found in the reviewed dictionary; they remain unassigned.", call. = FALSE)
  }
  x$data[[".nma_node_fine"]] <- dictionary[[fine_col]][key]
  x$data[[".nma_node_coarse"]] <- dictionary[[coarse_col]][key]
  x$mapping$node_fine <- ".nma_node_fine"
  x$mapping$node_coarse <- ".nma_node_coarse"

  if (x$format == "contrast") {
    label_col2 <- label_col2 %||% x$mapping$treatment2
    .nma_assert_cols(x$data, label_col2)
    key2 <- match(as.character(x$data[[label_col2]]), as.character(dictionary[[dictionary_label]]))
    unmatched2 <- is.na(key2) & !is.na(x$data[[label_col2]])
    if (any(unmatched2)) {
      warning(sum(unmatched2), " second-treatment labels were not found in the reviewed dictionary; they remain unassigned.", call. = FALSE)
    }
    x$data[[".nma_node_fine2"]] <- dictionary[[fine_col]][key2]
    x$data[[".nma_node_coarse2"]] <- dictionary[[coarse_col]][key2]
    x$mapping$node_fine2 <- ".nma_node_fine2"
    x$mapping$node_coarse2 <- ".nma_node_coarse2"
  }

  x$audit <- rbind(
    x$audit,
    .nma_audit_row(
      step = "node_policy",
      action = "apply_reviewed_dictionary",
      field = label_col,
      before = "Original labels preserved",
      after = paste0("fine=", fine_col, "; coarse=", coarse_col),
      rows = sum(!is.na(key)) + if (x$format == "contrast") sum(!is.na(key2)) else 0L,
      reason = reason %||% "Reviewed fine/coarse node dictionary",
      automatic = FALSE
    )
  )
  x
}

#' Split a network dataset into scientifically explicit strata
#'
#' This helper is useful when one input file contains multiple outcomes, time points,
#' populations, units, or other strata that should not be combined into a single NMA.
#' It does not decide which strata should be pooled; it makes the separation explicit.
#'
#' @param x An `nmaflow_data` object.
#' @param by Character vector of columns defining strata.
#' @param drop_empty Drop factor levels with no observations.
#' @return Named list of `nmaflow_data` objects.
#' @export
nma_split_strata <- function(x, by, drop_empty = TRUE) {
  x <- .nma_unwrap(x)
  .nma_assert_cols(x$data, by)
  key <- interaction(x$data[by], drop = drop_empty, sep = " | ")
  pieces <- split(seq_len(nrow(x$data)), key, drop = drop_empty)
  out <- lapply(names(pieces), function(nm) {
    z <- x
    z$data <- x$data[pieces[[nm]], , drop = FALSE]
    z$audit <- rbind(z$audit, .nma_audit_row(
      step = "stratify",
      action = "split_scientific_stratum",
      field = paste(by, collapse = ","),
      before = "Combined source dataset",
      after = nm,
      rows = nrow(z$data),
      reason = "Explicit scientific stratum",
      automatic = FALSE
    ))
    z
  })
  names(out) <- names(pieces)
  out
}

#' Compare network geometry before and after reviewed data operations
#'
#' @param before,after `nmaflow_data` objects.
#' @param node Node definition: `treatment`, `coarse`, or `fine`.
#' @return A diagnostic list comparing components, nodes, edges, and warnings.
#' @export
nma_network_repair_report <- function(before, after, node = c("treatment", "coarse", "fine")) {
  node <- match.arg(node)
  nb <- nma_network(before, node = node)
  na <- nma_network(after, node = node)
  vb <- nma_validate(before, strict = FALSE)
  va <- nma_validate(after, strict = FALSE)
  out <- list(
    node = node,
    before = data.frame(
      n_rows = nrow(.nma_unwrap(before)$data),
      n_nodes = nrow(nb$nodes),
      n_edges = nrow(nb$edges),
      n_components = length(nb$components),
      stringsAsFactors = FALSE
    ),
    after = data.frame(
      n_rows = nrow(.nma_unwrap(after)$data),
      n_nodes = nrow(na$nodes),
      n_edges = nrow(na$edges),
      n_components = length(na$components),
      stringsAsFactors = FALSE
    ),
    issues_before = vb$issues,
    issues_after = va$issues,
    audit = .nma_unwrap(after)$audit,
    quarantine = .nma_unwrap(after)$quarantine
  )
  class(out) <- c("nmaflow_repair_report", "list")
  out
}

#' @export
print.nmaflow_repair_report <- function(x, ...) {
  cat("nmaFlow network repair report\n")
  cat("Node level:", x$node, "\n\nBefore:\n")
  print(x$before, row.names = FALSE)
  cat("\nAfter:\n")
  print(x$after, row.names = FALSE)
  cat("\nAudit actions:", nrow(x$audit), "| Quarantined rows:", nrow(x$quarantine), "\n")
  invisible(x)
}
