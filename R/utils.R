`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

.nma_require <- function(pkg, feature = pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      sprintf("Feature '%s' requires optional package '%s'. Install it and retry.", feature, pkg),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

.nma_check_backend_result <- function(x, label = "backend") {
  engine <- if (is.list(x)) x$engine else NULL
  if (is.character(engine) && length(engine) == 1L && grepl("failed|error", engine, ignore.case = TRUE)) {
    stop(sprintf("%s failed: %s", label, engine), call. = FALSE)
  }
  if (inherits(x, "nma_nodesplit_df") && !is.null(x$model)) {
    valid <- vapply(x$model, function(z) {
      s <- z$stanfit
      inherits(s, "stanfit") && !is.null(s@sim$fnames_oi) && length(s@sim$fnames_oi) > 0L
    }, logical(1))
    if (any(!valid)) {
      stop(sprintf("%s failed: one or more node-split Stan models contain no posterior samples.", label), call. = FALSE)
    }
  }
  invisible(x)
}

.nma_assert_data <- function(data) {
  if (!is.data.frame(data)) stop("`data` must be a data.frame.", call. = FALSE)
  if (!nrow(data)) stop("`data` has zero rows.", call. = FALSE)
  invisible(TRUE)
}

.nma_assert_cols <- function(data, cols, label = "data") {
  cols <- unique(cols[!is.na(cols) & nzchar(cols)])
  miss <- setdiff(cols, names(data))
  if (length(miss)) {
    stop(sprintf("Missing column(s) in %s: %s", label, paste(miss, collapse = ", ")), call. = FALSE)
  }
  invisible(TRUE)
}

.nma_norm_text <- function(x) {
  x <- trimws(as.character(x))
  x <- gsub("[[:space:]]+", " ", x)
  x
}

.nma_norm_key <- function(x) {
  x <- tolower(.nma_norm_text(x))
  x <- iconv(x, to = "ASCII//TRANSLIT")
  x <- gsub("[^a-z0-9]+", " ", x)
  trimws(gsub("[[:space:]]+", " ", x))
}

.nma_empty_issues <- function() {
  data.frame(
    severity = character(),
    code = character(),
    study = character(),
    row = integer(),
    field = character(),
    message = character(),
    action = character(),
    stringsAsFactors = FALSE
  )
}

.nma_issue <- function(severity, code, message, study = NA_character_, row = NA_integer_,
                       field = NA_character_, action = "review") {
  data.frame(
    severity = severity,
    code = code,
    study = as.character(study),
    row = as.integer(row),
    field = field,
    message = message,
    action = action,
    stringsAsFactors = FALSE
  )
}

.nma_bind_issues <- function(...) {
  xs <- list(...)
  xs <- xs[vapply(xs, function(x) is.data.frame(x) && nrow(x) > 0L, logical(1))]
  if (!length(xs)) return(.nma_empty_issues())
  out <- do.call(rbind, xs)
  rownames(out) <- NULL
  out
}

.nma_new <- function(class, ..., call = match.call()) {
  structure(list(..., call = call), class = c(class, "nmaflow_object"))
}

.nma_audit_row <- function(step, action, field = NA_character_, before = NA_character_,
                           after = NA_character_, rows = NA_integer_, reason = NA_character_,
                           automatic = FALSE) {
  data.frame(
    timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE),
    step = step,
    action = action,
    field = field,
    before = before,
    after = after,
    rows = as.integer(rows),
    reason = reason,
    automatic = isTRUE(automatic),
    stringsAsFactors = FALSE
  )
}

.nma_audit_bind <- function(x, y) {
  if (is.null(x) || !nrow(x)) return(y)
  if (is.null(y) || !nrow(y)) return(x)
  out <- rbind(x, y)
  rownames(out) <- NULL
  out
}

.nma_study_arms <- function(data, study, treatment) {
  split(as.character(data[[treatment]]), as.character(data[[study]]))
}

.nma_edges_from_arms <- function(data, study, treatment) {
  arms <- .nma_study_arms(data, study, treatment)
  out <- lapply(names(arms), function(s) {
    trts <- unique(arms[[s]])
    trts <- trts[!is.na(trts) & nzchar(trts)]
    if (length(trts) < 2L) return(NULL)
    cmb <- utils::combn(sort(trts), 2L)
    data.frame(study = s, trt1 = cmb[1L, ], trt2 = cmb[2L, ], stringsAsFactors = FALSE)
  })
  out <- out[!vapply(out, is.null, logical(1))]
  if (!length(out)) return(data.frame(study = character(), trt1 = character(), trt2 = character()))
  do.call(rbind, out)
}

.nma_components <- function(nodes, edges) {
  nodes <- unique(as.character(nodes))
  nodes <- nodes[!is.na(nodes) & nzchar(nodes)]
  if (!length(nodes)) return(list())
  adj <- stats::setNames(vector("list", length(nodes)), nodes)
  if (nrow(edges)) {
    for (i in seq_len(nrow(edges))) {
      a <- as.character(edges$trt1[i]); b <- as.character(edges$trt2[i])
      adj[[a]] <- unique(c(adj[[a]], b)); adj[[b]] <- unique(c(adj[[b]], a))
    }
  }
  seen <- stats::setNames(rep(FALSE, length(nodes)), nodes)
  comps <- list()
  for (start in nodes) {
    if (seen[[start]]) next
    queue <- start; comp <- character()
    while (length(queue)) {
      v <- queue[[1L]]; queue <- queue[-1L]
      if (seen[[v]]) next
      seen[[v]] <- TRUE; comp <- c(comp, v)
      nxt <- adj[[v]]
      if (length(nxt)) queue <- c(queue, nxt[!seen[nxt]])
    }
    comps[[length(comps) + 1L]] <- sort(unique(comp))
  }
  comps
}

.nma_edge_summary <- function(edges) {
  if (!nrow(edges)) return(data.frame(trt1 = character(), trt2 = character(), studies = integer()))
  key <- paste(pmin(edges$trt1, edges$trt2), pmax(edges$trt1, edges$trt2), sep = "|||" )
  tab <- table(key)
  parts <- strsplit(names(tab), "\\|\\|\\|")
  data.frame(
    trt1 = vapply(parts, `[[`, character(1), 1L),
    trt2 = vapply(parts, `[[`, character(1), 2L),
    studies = as.integer(tab),
    stringsAsFactors = FALSE
  )
}

.nma_detect_format <- function(data, study, treatment, treatment2 = NULL) {
  if (!is.null(treatment2) && treatment2 %in% names(data)) return("contrast")
  if (study %in% names(data) && treatment %in% names(data)) return("arm")
  "unknown"
}
