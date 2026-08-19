#' List nmaFlow instructional blocks
#'
#' @param pattern Optional regular expression used to filter module or function names.
#' @return A data frame describing the package blocks.
#' @export
nma_blocks <- function(pattern = NULL) {
  x <- data.frame(
    block = 1:80,
    module = c(
      rep("Data architecture and network formation", 10),
      rep("Effects and network preflight", 6),
      rep("Frequentist NMA", 7),
      rep("Bayesian NMA", 7),
      rep("Bootstrap and uncertainty", 5),
      rep("Inconsistency and transitivity", 6),
      rep("Robustness and sensitivity", 5),
      rep("Meta-regression", 4),
      rep("Dose-response and nonlinear NMA", 5),
      rep("Complex evidence structures", 7),
      rep("Ranking and decision support", 4),
      rep("Publication output", 5),
      rep("Teaching, simulation and interoperability", 9)
    ),
    primary_function = c(
      "nma_data", "nma_validate", "nma_network", "nma_preflight", "nma_harmonize_candidates",
      "nma_apply_mapping", "nma_repair", "nma_reconnect_plan", "nma_quarantine", "nma_audit",
      "nma_effects", "nma_measure", "nma_direction", "nma_dispersion", "nma_multiarm_check", "nma_readiness",
      "nma_fit_freq", "nma_fit_nma", "nma_fit_additive", "nma_predict", "nma_direct", "nma_contribution", "nma_freq_summary",
      "nma_make_multinma", "nma_fit_bayes", "nma_prior", "nma_nodesplit_bayes", "nma_pp_check", "nma_mcmc_diagnose", "nma_bayes_summary",
      "nma_boot", "nma_boot_coef", "nma_boot_rank", "nma_boot_predict", "nma_boot_diagnose",
      "nma_transitivity", "nma_inconsistency", "nma_nodesplit", "nma_design_by_treatment", "nma_netheat", "nma_loop_check",
      "nma_influence", "nma_outliers", "nma_sensitivity", "nma_missing_sensitivity", "nma_zero_event_sensitivity",
      "nma_regression", "nma_effect_modifier", "nma_regression_plot", "nma_regression_compare",
      "nma_dose", "nma_nonlinear", "nma_dose_curve", "nma_dose_optimum", "nma_dose_compare",
      "nma_component", "nma_time", "nma_ipd", "nma_mlnmr", "nma_crossdesign", "nma_multivariate", "nma_sem",
      "nma_rank", "nma_rank_uncertainty", "nma_utility", "nma_threshold",
      "nma_plot", "nma_table", "nma_league", "nma_sof", "nma_report",
      "nma_simulate", "nma_data_example", "nma_tour", "nma_explain", "nma_capabilities",
      "nma_export", "nma_reproduce", "nma_session", "nma_workflow"
    ),
    stringsAsFactors = FALSE
  )
  if (!is.null(pattern)) {
    keep <- grepl(pattern, x$module, ignore.case = TRUE) | grepl(pattern, x$primary_function, ignore.case = TRUE)
    x <- x[keep, , drop = FALSE]
  }
  rownames(x) <- NULL
  x
}

#' Report optional backend availability
#'
#' @return A data frame with package availability and intended role.
#' @export
nma_capabilities <- function() {
  pkgs <- c("netmeta", "NMA", "multinma", "gemtc", "nmaINLA", "MBNMAdose", "netdose",
            "MBNMAtime", "crossnma", "mixmeta", "metaSEM", "lavaan", "ggplot2", "igraph",
            "gt", "flextable", "reticulate")
  roles <- c(
    "frequentist NMA", "multivariate/frequentist alternative", "Bayesian NMA and ML-NMR",
    "Bayesian JAGS compatibility", "Bayesian INLA", "Bayesian dose-response", "frequentist dose-response",
    "time-course NMA", "cross-design/cross-format", "multivariate meta-analysis", "MASEM/SEM",
    "second-stage SEM", "publication graphics", "network geometry graphics", "publication tables",
    "Word-ready tables", "optional Python bridge"
  )
  data.frame(
    package = pkgs,
    available = vapply(pkgs, requireNamespace, logical(1), quietly = TRUE),
    role = roles,
    stringsAsFactors = FALSE
  )
}
