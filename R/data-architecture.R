#' Create an auditable nmaFlow data object
#'
#' @param data A data frame in arm-level or contrast-level format.
#' @param study Name of the study identifier column.
#' @param treatment Name of the treatment column for arm-level data or first treatment for contrast data.
#' @param treatment2 Optional second treatment column for contrast-level data.
#' @param arm Optional arm identifier column.
#' @param outcome Optional outcome identifier column.
#' @param time Optional follow-up/timepoint column.
#' @param unit Optional outcome unit column.
#' @param node_coarse Optional primary-analysis treatment node column.
#' @param node_fine Optional fine treatment node column for dose/component sensitivity analyses.
#' @param node_coarse2,node_fine2 Optional second-treatment coarse/fine node columns for contrast-level data.
#' @param population Optional analysis-population column (for example ITT or PP).
#' @param format One of `auto`, `arm`, or `contrast`.
#' @param metadata Optional named list stored with the object.
#' @return An object of class `nmaflow_data`.
#' @export
nma_data <- function(data, study = "study", treatment = "treatment", treatment2 = NULL,
                     arm = NULL, outcome = NULL, time = NULL, unit = NULL,
                     node_coarse = NULL, node_fine = NULL, node_coarse2 = NULL, node_fine2 = NULL, population = NULL,
                     format = c("auto", "arm", "contrast"), metadata = list()) {
  .nma_assert_data(data)
  format <- match.arg(format)
  .nma_assert_cols(data, c(study, treatment, treatment2, arm, outcome, time, unit,
                           node_coarse, node_fine, node_coarse2, node_fine2, population))
  if (format == "auto") format <- .nma_detect_format(data, study, treatment, treatment2)
  if (format == "unknown") stop("Could not infer data format. Supply `format` explicitly.", call. = FALSE)

  mapping <- list(
    study = study, treatment = treatment, treatment2 = treatment2, arm = arm,
    outcome = outcome, time = time, unit = unit, node_coarse = node_coarse,
    node_fine = node_fine, node_coarse2 = node_coarse2, node_fine2 = node_fine2, population = population
  )
  audit <- .nma_audit_row("import", "create_nma_data", rows = nrow(data),
                          reason = sprintf("Imported %s-level evidence", format), automatic = FALSE)
  .nma_new(
    "nmaflow_data",
    data = data,
    raw_data = data,
    quarantine = data[0, , drop = FALSE],
    audit = audit,
    format = format,
    mapping = mapping,
    metadata = metadata
  )
}

#' @export
print.nmaflow_data <- function(x, ...) {
  cat("nmaFlow data object\n")
  cat("  format:", x$format, "\n")
  cat("  rows:", nrow(x$data), "\n")
  cat("  quarantined:", nrow(x$quarantine), "\n")
  st <- x$mapping$study; tr <- x$mapping$treatment
  if (!is.null(st) && st %in% names(x$data)) cat("  studies:", length(unique(x$data[[st]])), "\n")
  if (!is.null(tr) && tr %in% names(x$data)) cat("  treatments:", length(unique(x$data[[tr]])), "\n")
  invisible(x)
}

.nma_unwrap <- function(x) {
  if (inherits(x, "nmaflow_data")) return(x)
  stop("Expected an `nmaflow_data` object created by nma_data().", call. = FALSE)
}

#' Validate data architecture before network meta-analysis
#'
#' The validator identifies structural errors without modifying the data. It checks
#' identifiers, duplicated rows/arms, one-arm studies, duplicated treatment nodes within
#' a study, incompatible outcome/time/unit strata, missing dispersion, and network geometry.
#'
#' @param x An `nmaflow_data` object.
#' @param n Optional sample-size column.
#' @param mean Optional mean column.
#' @param sd Optional standard-deviation column.
#' @param se Optional standard-error column.
#' @param event Optional event-count column.
#' @param total Optional denominator column for binary data.
#' @param strict If `TRUE`, readiness is false when any error is present.
#' @return An object of class `nmaflow_validation`.
#' @export
nma_validate <- function(x, n = NULL, mean = NULL, sd = NULL, se = NULL,
                         event = NULL, total = NULL, strict = TRUE) {
  x <- .nma_unwrap(x)
  d <- x$data; m <- x$mapping
  .nma_assert_cols(d, c(n, mean, sd, se, event, total))
  issues <- .nma_empty_issues()
  st <- m$study; tr <- m$treatment

  bad_study <- which(is.na(d[[st]]) | !nzchar(trimws(as.character(d[[st]]))))
  if (length(bad_study)) issues <- .nma_bind_issues(issues, do.call(rbind, lapply(bad_study, function(i)
    .nma_issue("error", "missing_study_id", "Study identifier is missing.", row = i, field = st, action = "quarantine"))))

  bad_trt <- which(is.na(d[[tr]]) | !nzchar(trimws(as.character(d[[tr]]))))
  if (length(bad_trt)) issues <- .nma_bind_issues(issues, do.call(rbind, lapply(bad_trt, function(i)
    .nma_issue("error", "missing_treatment", "Treatment node is missing.", study = d[[st]][i], row = i, field = tr, action = "review"))))

  dup_rows <- which(duplicated(d))
  if (length(dup_rows)) issues <- .nma_bind_issues(issues, do.call(rbind, lapply(dup_rows, function(i)
    .nma_issue("warning", "duplicate_row", "Exact duplicate row detected.", study = d[[st]][i], row = i, action = "quarantine_duplicate"))))

  if (x$format == "arm") {
    strata_cols <- c(m$outcome, m$time, m$population)
    strata_cols <- strata_cols[!vapply(strata_cols, is.null, logical(1))]
    key_cols <- c(st, strata_cols, if (!is.null(m$arm)) m$arm else tr)
    key <- do.call(paste, c(d[key_cols], sep = "|||"))
    dups <- which(duplicated(key) | duplicated(key, fromLast = TRUE))
    if (length(dups)) {
      for (i in dups) {
        issues <- .nma_bind_issues(issues, .nma_issue(
          "warning", "duplicate_arm_or_node",
          "The same study/stratum arm or treatment node occurs more than once. Review whether rows are duplicated outcomes, repeated measurements, or a coding error.",
          study = d[[st]][i], row = i, field = if (!is.null(m$arm)) m$arm else tr, action = "review"
        ))
      }
    }

    arms <- split(d, as.character(d[[st]]))
    for (s in names(arms)) {
      z <- arms[[s]]
      nts <- length(unique(as.character(z[[tr]])))
      if (nts < 2L) issues <- .nma_bind_issues(issues, .nma_issue(
        "error", "single_arm_study", "Study contributes fewer than two distinct treatment nodes.", study = s,
        field = tr, action = "exclude_from_contrast_nma_or_review"
      ))
      if (any(duplicated(as.character(z[[tr]])))) {
        issues <- .nma_bind_issues(issues, .nma_issue(
          "warning", "same_node_within_study",
          "A study contains repeated rows assigned to the same treatment node. If arms are scientifically distinct, refine node coding; if they are repeated outcome/time rows, stratify before network construction.",
          study = s, field = tr, action = "review_node_policy"
        ))
      }
    }
  }

  if (!is.null(m$outcome)) {
    bad <- which(is.na(d[[m$outcome]]) | !nzchar(trimws(as.character(d[[m$outcome]]))))
    if (length(bad)) issues <- .nma_bind_issues(issues, .nma_issue(
      "error", "missing_outcome", sprintf("%d row(s) have missing outcome identifiers.", length(bad)), field = m$outcome, action = "review"
    ))
  }

  if (!is.null(m$unit) && !is.null(m$outcome)) {
    ou <- split(as.character(d[[m$unit]]), as.character(d[[m$outcome]]))
    for (o in names(ou)) {
      u <- unique(ou[[o]][!is.na(ou[[o]]) & nzchar(ou[[o]])])
      if (length(u) > 1L) issues <- .nma_bind_issues(issues, .nma_issue(
        "warning", "mixed_units", sprintf("Outcome '%s' has multiple units: %s", o, paste(u, collapse = ", ")),
        field = m$unit, action = "convert_with_explicit_rule_or_stratify"
      ))
    }
  }

  if (!is.null(m$time) && !is.null(m$outcome)) {
    ot <- split(as.character(d[[m$time]]), as.character(d[[m$outcome]]))
    for (o in names(ot)) {
      tt <- unique(ot[[o]][!is.na(ot[[o]]) & nzchar(ot[[o]])])
      if (length(tt) > 1L) issues <- .nma_bind_issues(issues, .nma_issue(
        "info", "multiple_timepoints", sprintf("Outcome '%s' contains %d timepoints. Build separate networks or use a time-course model.", o, length(tt)),
        field = m$time, action = "stratify_or_time_model"
      ))
    }
  }

  if (!is.null(n)) {
    bad <- which(is.na(d[[n]]) | !is.finite(d[[n]]) | d[[n]] <= 0)
    if (length(bad)) issues <- .nma_bind_issues(issues, .nma_issue(
      "error", "invalid_sample_size", sprintf("%d row(s) have missing/non-positive sample size.", length(bad)), field = n, action = "review"
    ))
  }
  if (!is.null(sd)) {
    bad <- which(is.na(d[[sd]]) | !is.finite(d[[sd]]) | d[[sd]] <= 0)
    if (length(bad)) issues <- .nma_bind_issues(issues, .nma_issue(
      "warning", "missing_or_invalid_sd", sprintf("%d row(s) lack a positive SD. Use documented recovery rules or sensitivity analysis.", length(bad)), field = sd, action = "recover_or_mark_unavailable"
    ))
  }
  if (!is.null(se)) {
    bad <- which(!is.na(d[[se]]) & (!is.finite(d[[se]]) | d[[se]] <= 0))
    if (length(bad)) issues <- .nma_bind_issues(issues, .nma_issue(
      "error", "invalid_se", sprintf("%d row(s) have non-positive SE.", length(bad)), field = se, action = "review"
    ))
  }
  if (!is.null(event) && !is.null(total)) {
    bad <- which(is.na(d[[event]]) | is.na(d[[total]]) | d[[event]] < 0 | d[[total]] <= 0 | d[[event]] > d[[total]])
    if (length(bad)) issues <- .nma_bind_issues(issues, .nma_issue(
      "error", "invalid_binomial_counts", sprintf("%d row(s) have invalid event/total counts.", length(bad)), action = "review"
    ))
  }

  geom <- tryCatch(nma_network(x), error = function(e) NULL)
  if (!is.null(geom) && length(geom$components) > 1L) {
    issues <- .nma_bind_issues(issues, .nma_issue(
      "error", "disconnected_network",
      sprintf("Network has %d disconnected components. Do not force-connect them; review node coding, stratification and evidence eligibility.", length(geom$components)),
      action = "run_nma_reconnect_plan"
    ))
  }

  sev <- table(factor(issues$severity, levels = c("error", "warning", "info")))
  ready <- if (strict) unname(sev[["error"]] == 0L) else TRUE
  .nma_new("nmaflow_validation", issues = issues, counts = sev, ready = ready, network = geom, mapping = m)
}

#' @export
print.nmaflow_validation <- function(x, ...) {
  cat("nmaFlow data validation\n")
  cat("  ready:", if (isTRUE(x$ready)) "YES" else "NO", "\n")
  cat("  errors:", unname(x$counts[["error"]]), " warnings:", unname(x$counts[["warning"]]),
      " info:", unname(x$counts[["info"]]), "\n")
  if (nrow(x$issues)) print(utils::head(x$issues, 12L), row.names = FALSE)
  invisible(x)
}

#' Build and summarize network geometry
#'
#' @param x An `nmaflow_data` object.
#' @param node Which treatment representation to use: `treatment`, `coarse`, or `fine`.
#' @param subset Optional logical vector used to construct a stratum-specific network.
#' @return An object of class `nmaflow_network`.
#' @export
nma_network <- function(x, node = c("treatment", "coarse", "fine"), subset = NULL) {
  x <- .nma_unwrap(x); node <- match.arg(node)
  d <- x$data; m <- x$mapping
  if (!is.null(subset)) d <- d[subset, , drop = FALSE]
  tr <- switch(node,
    treatment = m$treatment,
    coarse = m$node_coarse %||% stop("No `node_coarse` column was declared.", call. = FALSE),
    fine = m$node_fine %||% stop("No `node_fine` column was declared.", call. = FALSE)
  )
  .nma_assert_cols(d, c(m$study, tr))
  if (x$format == "contrast") {
    tr2 <- switch(node,
      treatment = m$treatment2,
      coarse = m$node_coarse2 %||% stop("Contrast-level coarse networks require `node_coarse2`.", call. = FALSE),
      fine = m$node_fine2 %||% stop("Contrast-level fine networks require `node_fine2`.", call. = FALSE)
    )
    .nma_assert_cols(d, tr2)
    edges <- data.frame(study = as.character(d[[m$study]]), trt1 = as.character(d[[tr]]),
                        trt2 = as.character(d[[tr2]]), stringsAsFactors = FALSE)
  } else {
    edges <- .nma_edges_from_arms(d, m$study, tr)
  }
  edges <- edges[!is.na(edges$trt1) & !is.na(edges$trt2) & edges$trt1 != edges$trt2, , drop = FALSE]
  nodes <- sort(unique(c(edges$trt1, edges$trt2, as.character(d[[tr]]))))
  nodes <- nodes[!is.na(nodes) & nzchar(nodes)]
  comps <- .nma_components(nodes, edges)
  es <- .nma_edge_summary(edges)
  deg <- stats::setNames(integer(length(nodes)), nodes)
  if (nrow(es)) {
    for (i in seq_len(nrow(es))) {
      deg[es$trt1[i]] <- deg[es$trt1[i]] + es$studies[i]
      deg[es$trt2[i]] <- deg[es$trt2[i]] + es$studies[i]
    }
  }
  .nma_new(
    "nmaflow_network",
    node_type = node,
    nodes = data.frame(treatment = nodes, weighted_degree = as.integer(deg[nodes]), stringsAsFactors = FALSE),
    edges = es,
    study_edges = edges,
    components = comps,
    connected = length(comps) <= 1L,
    n_studies = length(unique(as.character(d[[m$study]])))
  )
}

#' @export
print.nmaflow_network <- function(x, ...) {
  cat("nmaFlow network geometry\n")
  cat("  treatments:", nrow(x$nodes), "\n")
  cat("  studies:", x$n_studies, "\n")
  cat("  observed comparisons:", nrow(x$edges), "\n")
  cat("  connected:", if (x$connected) "YES" else "NO", "\n")
  if (!x$connected) {
    cat("  components:\n")
    for (i in seq_along(x$components)) cat("   ", i, ":", paste(x$components[[i]], collapse = ", "), "\n")
  }
  invisible(x)
}

#' Preflight check for NMA readiness
#'
#' @param x An `nmaflow_data` object.
#' @param ... Passed to [nma_validate()].
#' @return An object of class `nmaflow_preflight` with GREEN/YELLOW/RED status.
#' @export
nma_preflight <- function(x, ...) {
  val <- nma_validate(x, ...)
  status <- if (unname(val$counts[["error"]]) > 0L) "RED" else if (unname(val$counts[["warning"]]) > 0L) "YELLOW" else "GREEN"
  .nma_new("nmaflow_preflight", status = status, validation = val, network = val$network)
}

#' @export
print.nmaflow_preflight <- function(x, ...) {
  cat("nmaFlow network preflight:", x$status, "\n")
  print(x$validation)
  invisible(x)
}

#' Generate treatment harmonization candidates without modifying data
#'
#' Candidate generation uses conservative normalized-string similarity. It is a review aid,
#' not an automatic node-merging procedure.
#'
#' @param x An `nmaflow_data` object.
#' @param max_distance Maximum normalized edit distance (0 to 1).
#' @param node Treatment representation to inspect.
#' @return A data frame of candidate mappings requiring review.
#' @export
nma_harmonize_candidates <- function(x, max_distance = 0.18,
                                     node = c("treatment", "coarse", "fine")) {
  x <- .nma_unwrap(x); node <- match.arg(node); m <- x$mapping
  col <- switch(node, treatment = m$treatment, coarse = m$node_coarse, fine = m$node_fine)
  if (is.null(col)) stop(sprintf("No %s node column declared.", node), call. = FALSE)
  vals <- as.character(x$data[[col]])
  if (x$format == "contrast") {
    col2 <- switch(node, treatment = m$treatment2, coarse = m$node_coarse2, fine = m$node_fine2)
    if (!is.null(col2) && col2 %in% names(x$data)) vals <- c(vals, as.character(x$data[[col2]]))
  }
  vals <- sort(unique(vals))
  vals <- vals[!is.na(vals) & nzchar(vals)]
  if (length(vals) < 2L) return(data.frame())
  cmb <- utils::combn(vals, 2L)
  out <- vector("list", ncol(cmb))
  for (i in seq_len(ncol(cmb))) {
    a <- cmb[1L, i]; b <- cmb[2L, i]
    ak <- .nma_norm_key(a); bk <- .nma_norm_key(b)
    den <- max(nchar(ak), nchar(bk), 1L)
    dist <- utils::adist(ak, bk)[1L] / den
    exact_key <- identical(ak, bk)
    if (exact_key || dist <= max_distance) {
      out[[i]] <- data.frame(from = b, to_candidate = a, normalized_distance = dist,
                             exact_normalized_key = exact_key,
                             recommendation = if (exact_key) "safe_formatting_candidate" else "manual_review_only",
                             stringsAsFactors = FALSE)
    }
  }
  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) return(data.frame())
  ans <- do.call(rbind, out); rownames(ans) <- NULL; ans[order(ans$normalized_distance), , drop = FALSE]
}

#' Apply an explicitly reviewed node or category mapping
#'
#' @param x An `nmaflow_data` object.
#' @param mapping A two-column data frame or named character vector mapping old to reviewed values.
#' @param field Field to modify (`treatment`, `coarse`, `fine`, `outcome`, `time`, `unit`, or a literal column name).
#' @param reviewed Must be `TRUE`; this guard prevents accidental semantic auto-merging.
#' @param reason Required free-text rationale recorded in the audit log.
#' @return Updated `nmaflow_data` object.
#' @export
nma_apply_mapping <- function(x, mapping, field = "treatment", reviewed = FALSE, reason = NULL) {
  x <- .nma_unwrap(x)
  if (!isTRUE(reviewed)) stop("Set `reviewed = TRUE` only after scientific review of the mapping.", call. = FALSE)
  if (is.null(reason) || !nzchar(reason)) stop("A non-empty `reason` is required for reviewed mappings.", call. = FALSE)
  col <- switch(field,
    treatment = x$mapping$treatment,
    coarse = x$mapping$node_coarse,
    fine = x$mapping$node_fine,
    coarse2 = x$mapping$node_coarse2,
    fine2 = x$mapping$node_fine2,
    outcome = x$mapping$outcome,
    time = x$mapping$time,
    unit = x$mapping$unit,
    field
  )
  if (is.null(col) || !col %in% names(x$data)) stop("Requested mapping field is not available.", call. = FALSE)
  if (is.data.frame(mapping)) {
    if (ncol(mapping) < 2L) stop("`mapping` data frame needs at least two columns: from and to.", call. = FALSE)
    mp <- stats::setNames(as.character(mapping[[2L]]), as.character(mapping[[1L]]))
  } else {
    mp <- mapping
    if (is.null(names(mp))) stop("Character mapping must be named: old = new.", call. = FALSE)
  }
  old <- as.character(x$data[[col]])
  hit <- !is.na(old) & old %in% names(mp)
  new <- old
  new[hit] <- unname(mp[old[hit]])
  x$data[[col]] <- new
  x$audit <- .nma_audit_bind(x$audit, .nma_audit_row(
    "mapping", "apply_reviewed_mapping", field = col,
    before = paste(unique(old[hit]), collapse = " | "),
    after = paste(unique(new[hit]), collapse = " | "), rows = sum(hit), reason = reason, automatic = FALSE
  ))
  x
}

#' Repair unambiguous data-architecture problems and quarantine unsafe rows
#'
#' @param x An `nmaflow_data` object.
#' @param trim_whitespace Normalize leading/trailing/repeated whitespace in declared character fields.
#' @param normalize_equivalent_spelling Merge labels only when lower-case transliterated punctuation-free keys are exactly identical.
#' @param quarantine_exact_duplicates Move exact duplicate rows after the first copy to `x$quarantine`.
#' @param explicit_mapping Optional reviewed mapping passed to [nma_apply_mapping()].
#' @param mapping_field Field for `explicit_mapping`.
#' @param mapping_reason Rationale for the explicit mapping.
#' @return Updated `nmaflow_data` object with audit log.
#' @export
nma_repair <- function(x, trim_whitespace = TRUE, normalize_equivalent_spelling = TRUE,
                       quarantine_exact_duplicates = TRUE, explicit_mapping = NULL,
                       mapping_field = "treatment", mapping_reason = NULL) {
  x <- .nma_unwrap(x)
  d <- x$data
  declared <- unique(unlist(x$mapping, use.names = FALSE))
  declared <- declared[!is.na(declared) & declared %in% names(d)]

  if (isTRUE(trim_whitespace)) {
    for (col in declared) {
      if (is.character(d[[col]]) || is.factor(d[[col]])) {
        old <- as.character(d[[col]]); new <- .nma_norm_text(old)
        changed <- !is.na(old) & old != new
        if (any(changed)) {
          d[[col]] <- new
          x$audit <- .nma_audit_bind(x$audit, .nma_audit_row(
            "repair", "normalize_whitespace", field = col, rows = sum(changed),
            reason = "Whitespace-only normalization; no semantic mapping", automatic = TRUE
          ))
        }
      }
    }
  }

  if (isTRUE(normalize_equivalent_spelling)) {
    for (field in c(x$mapping$treatment, x$mapping$treatment2, x$mapping$node_coarse, x$mapping$node_fine,
                    x$mapping$node_coarse2, x$mapping$node_fine2)) {
      if (is.null(field) || !field %in% names(d)) next
      vals <- unique(as.character(d[[field]])); vals <- vals[!is.na(vals)]
      keys <- .nma_norm_key(vals)
      groups <- split(vals, keys)
      groups <- groups[vapply(groups, length, integer(1)) > 1L]
      for (g in groups) {
        canonical <- g[[1L]]
        hit <- as.character(d[[field]]) %in% g
        before <- paste(g, collapse = " | ")
        d[[field]][hit] <- canonical
        x$audit <- .nma_audit_bind(x$audit, .nma_audit_row(
          "repair", "merge_formatting_equivalents", field = field, before = before, after = canonical,
          rows = sum(hit), reason = "Labels differ only by case, accents, punctuation or whitespace after strict normalization",
          automatic = TRUE
        ))
      }
    }
  }

  if (isTRUE(quarantine_exact_duplicates)) {
    dup <- duplicated(d)
    if (any(dup)) {
      x$quarantine <- rbind(x$quarantine, d[dup, , drop = FALSE])
      d <- d[!dup, , drop = FALSE]
      x$audit <- .nma_audit_bind(x$audit, .nma_audit_row(
        "repair", "quarantine_exact_duplicates", rows = sum(dup),
        reason = "Exact duplicate rows retained in quarantine rather than silently deleted", automatic = TRUE
      ))
    }
  }

  x$data <- d
  if (!is.null(explicit_mapping)) {
    x <- nma_apply_mapping(x, explicit_mapping, field = mapping_field, reviewed = TRUE,
                           reason = mapping_reason %||% "User-supplied reviewed mapping")
  }
  x
}

#' Diagnose disconnected components and propose review targets
#'
#' This function never edits treatment nodes or invents bridging comparisons. It searches for
#' suspiciously similar labels across disconnected components and reports study/outcome/time
#' strata that may explain fragmentation.
#'
#' @param x An `nmaflow_data` object.
#' @param max_distance Maximum normalized string distance used only for candidate review.
#' @param node Node representation used to diagnose the network.
#' @return An object of class `nmaflow_reconnect_plan`.
#' @export
nma_reconnect_plan <- function(x, max_distance = 0.25, node = c("treatment", "coarse", "fine")) {
  x <- .nma_unwrap(x); node <- match.arg(node)
  net <- nma_network(x, node = node)
  if (net$connected) return(.nma_new("nmaflow_reconnect_plan", connected = TRUE, components = net$components,
                                    candidates = data.frame(), recommendations = "Network is already connected."))
  comps <- net$components
  cand <- list(); k <- 0L
  if (length(comps) > 1L) {
    for (i in seq_len(length(comps) - 1L)) for (j in (i + 1L):length(comps)) {
      for (a in comps[[i]]) for (b in comps[[j]]) {
        ak <- .nma_norm_key(a); bk <- .nma_norm_key(b)
        dist <- utils::adist(ak, bk)[1L] / max(nchar(ak), nchar(bk), 1L)
        if (dist <= max_distance) {
          k <- k + 1L
          cand[[k]] <- data.frame(component1 = i, node1 = a, component2 = j, node2 = b,
                                  normalized_distance = dist,
                                  interpretation = "possible coding/alias issue; review source definitions before any merge",
                                  stringsAsFactors = FALSE)
        }
      }
    }
  }
  candidates <- if (length(cand)) do.call(rbind, cand) else data.frame()
  rec <- c(
    "Check whether near-duplicate node labels represent the same intervention; apply only reviewed mappings.",
    "Check outcome, timepoint, population and unit stratification before combining rows.",
    "Check whether active co-interventions or dose distinctions should remain separate scientific nodes.",
    "Check eligibility/extraction for omitted comparator arms.",
    "If components remain scientifically distinct after review, analyze them separately; do not fabricate a bridge."
  )
  .nma_new("nmaflow_reconnect_plan", connected = FALSE, components = comps, candidates = candidates,
           recommendations = rec)
}

#' @export
print.nmaflow_reconnect_plan <- function(x, ...) {
  cat("nmaFlow reconnection diagnostic\n")
  cat("  connected:", if (x$connected) "YES" else "NO", "\n")
  if (!x$connected) {
    cat("  candidate coding issues:", nrow(x$candidates), "\n")
    if (nrow(x$candidates)) print(utils::head(x$candidates, 12L), row.names = FALSE)
    cat("  recommendations:\n -", paste(x$recommendations, collapse = "\n - "), "\n")
  }
  invisible(x)
}

#' Quarantine selected rows with a documented reason
#'
#' @param x An `nmaflow_data` object.
#' @param rows Integer or logical rows to quarantine.
#' @param reason Required reason.
#' @return Updated data object.
#' @export
nma_quarantine <- function(x, rows, reason) {
  x <- .nma_unwrap(x)
  if (missing(reason) || !nzchar(reason)) stop("A quarantine reason is required.", call. = FALSE)
  idx <- seq_len(nrow(x$data))[rows]
  idx <- unique(idx[is.finite(idx) & idx >= 1L & idx <= nrow(x$data)])
  if (!length(idx)) return(x)
  x$quarantine <- rbind(x$quarantine, x$data[idx, , drop = FALSE])
  x$data <- x$data[-idx, , drop = FALSE]
  x$audit <- .nma_audit_bind(x$audit, .nma_audit_row("quarantine", "quarantine_rows", rows = length(idx),
                                                      reason = reason, automatic = FALSE))
  x
}

#' Retrieve the audit trail
#'
#' @param x An `nmaflow_data` object.
#' @return Audit data frame.
#' @export
nma_audit <- function(x) {
  x <- .nma_unwrap(x)
  x$audit
}
