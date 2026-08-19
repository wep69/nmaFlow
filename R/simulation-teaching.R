#' Simulate agronomic network meta-analysis data
#'
#' @param domain Agronomic teaching domain.
#' @param studies Number of studies.
#' @param seed Random seed.
#' @param outcome_type `continuous` or `binary`.
#' @param multiarm_fraction Approximate probability of a three-arm study.
#' @param heterogeneity Between-study SD on the outcome/link scale.
#' @param inconsistency Optional deterministic inconsistency shift added to one comparison design.
#' @return Arm-level data frame containing truth columns for method validation.
#' @export
nma_simulate <- function(domain = c("nitrogen", "irrigation", "biostimulants", "crop_protection",
                                    "soil_amendments", "multienvironment", "timecourse", "multioutcome", "inconsistent"),
                         studies = 30, seed = getOption("nmaFlow.seed", 260819),
                         outcome_type = c("continuous", "binary"), multiarm_fraction = 0.20,
                         heterogeneity = 0.12, inconsistency = 0) {
  inconsistency_missing <- missing(inconsistency)
  domain <- match.arg(domain); outcome_type <- match.arg(outcome_type)
  if (domain == "inconsistent" && inconsistency_missing) inconsistency <- 0.35
  if (studies < 3L) stop("Use at least three studies for an NMA teaching dataset.", call. = FALSE)
  set.seed(seed)
  defs <- switch(domain,
    nitrogen = list(trt = c("Control", "Urea", "Ammonium sulfate", "Urea+inhibitor", "Organic N"),
                    eff = c(0, 0.45, 0.38, 0.55, 0.30)),
    irrigation = list(trt = c("Rainfed", "Deficit 50%", "Deficit 75%", "Full irrigation", "Sensor-guided"),
                      eff = c(0, 0.25, 0.42, 0.55, 0.60)),
    biostimulants = list(trt = c("Control", "Seaweed", "Amino acids", "Humic acids", "Microbial", "Seaweed+Amino acids"),
                        eff = c(0, 0.22, 0.18, 0.15, 0.25, 0.35)),
    crop_protection = list(trt = c("Untreated", "Chemical A", "Chemical B", "Biological", "IPM"),
                           eff = c(0, -0.55, -0.45, -0.30, -0.62)),
    soil_amendments = list(trt = c("Control", "Lime", "Gypsum", "Biochar", "Lime+Gypsum", "Biochar+Lime"),
                           eff = c(0, 0.28, 0.18, 0.22, 0.40, 0.45)),
    multienvironment = list(trt = paste0("Genotype ", LETTERS[1:6]), eff = c(0, .16, .22, .10, .28, .19)),
    timecourse = list(trt = c("Control", "Treatment A", "Treatment B", "Integrated"), eff = c(0, .22, .30, .42)),
    multioutcome = list(trt = c("Control", "Practice A", "Practice B", "Integrated"), eff = c(0, .20, .28, .40)),
    inconsistent = list(trt = c("Control", "A", "B", "C"), eff = c(0, .25, .40, .55))
  )
  trts <- defs$trt; te <- stats::setNames(defs$eff, trts)
  rows <- list(); k <- 0L
  for (s in seq_len(studies)) {
    narm <- if (stats::runif(1) < multiarm_fraction) 3L else 2L
    # Anchor approximately half of studies; otherwise create head-to-head evidence.
    if (stats::runif(1) < 0.55) arms <- unique(c(trts[1L], sample(trts[-1L], narm - 1L)))
    else arms <- sample(trts, narm, replace = FALSE)
    if (length(arms) < 2L) next
    baseline <- stats::rnorm(1, if (domain == "crop_protection") 0 else 6.5, 0.6)
    u <- stats::rnorm(1, 0, heterogeneity)
    soil_pH <- stats::rnorm(1, 6.1, 0.55)
    rainfall <- max(100, stats::rnorm(1, 720, 180))
    year <- sample(2010:2026, 1)
    for (a in seq_along(arms)) {
      tr <- arms[a]; nn <- sample(20:80, 1)
      design_shift <- 0
      if (domain == "inconsistent" && all(c("A", "C") %in% arms) && tr == "C") design_shift <- inconsistency
      lp <- baseline + te[[tr]] + u + design_shift + 0.03 * (soil_pH - 6) + 0.0002 * (rainfall - 700)
      k <- k + 1L
      if (outcome_type == "continuous") {
        sd <- stats::runif(1, 0.65, 1.35)
        mean_obs <- stats::rnorm(1, lp, sd / sqrt(nn))
        rows[[k]] <- data.frame(study = sprintf("S%03d", s), arm = a, treatment = tr,
                                mean = mean_obs, sd = sd, n = nn, soil_pH = soil_pH,
                                rainfall = rainfall, year = year, true_effect = te[[tr]] + design_shift,
                                domain = domain, stringsAsFactors = FALSE)
      } else {
        p <- stats::plogis(lp - 6.5)
        ev <- stats::rbinom(1, nn, p)
        rows[[k]] <- data.frame(study = sprintf("S%03d", s), arm = a, treatment = tr,
                                event = ev, n = nn, soil_pH = soil_pH, rainfall = rainfall,
                                year = year, true_logodds_effect = te[[tr]] + design_shift,
                                domain = domain, stringsAsFactors = FALSE)
      }
    }
  }
  out <- do.call(rbind, rows)
  if (domain == "nitrogen") {
    dose_map <- c("Control" = 0, "Urea" = 100, "Ammonium sulfate" = 100,
                  "Urea+inhibitor" = 100, "Organic N" = 100)
    out$dose <- unname(dose_map[out$treatment]); out$agent <- out$treatment
    out$node_coarse <- sub("\\+inhibitor", "", out$treatment)
    out$node_fine <- paste0(out$treatment, "|", out$dose)
  }
  if (domain == "biostimulants") {
    out$node_coarse <- out$treatment
    out$components <- gsub("\\+", " + ", out$treatment)
  }
  if (domain == "soil_amendments") out$components <- gsub("\\+", " + ", out$treatment)
  if (domain == "timecourse") {
    out <- do.call(rbind, lapply(seq_len(nrow(out)), function(i) {
      times <- c(30, 60, 90, 120)
      z <- out[rep(i, length(times)), , drop = FALSE]
      z$time <- times
      z$mean <- z$mean + 0.004 * times + stats::rnorm(length(times), 0, z$sd / sqrt(z$n))
      z
    }))
  }
  if (domain == "multioutcome" && outcome_type == "continuous") {
    out$outcome <- "yield"; out$unit <- "t ha-1"
    q <- out; q$outcome <- "protein"; q$unit <- "%"; q$mean <- 10 + 0.5 * (q$mean - mean(q$mean)) + stats::rnorm(nrow(q), 0, .2)
    w <- out; w$outcome <- "WUE"; w$unit <- "kg m-3"; w$mean <- 1.5 + 0.15 * (w$mean - mean(w$mean)) + stats::rnorm(nrow(w), 0, .08)
    out <- rbind(out, q, w)
  }
  rownames(out) <- NULL
  attr(out, "simulation_truth") <- list(seed = seed, treatment_effects = te, heterogeneity = heterogeneity,
                                         inconsistency = inconsistency, domain = domain)
  out
}

#' Retrieve a frozen-style teaching dataset
#'
#' @param name One of the documented agronomic examples.
#' @return A reproducibly simulated data frame.
#' @export
nma_data_example <- function(name = c("agro_nitrogen", "agro_irrigation", "agro_biostimulants",
                                      "agro_crop_protection", "agro_soil_amendments", "agro_multienv",
                                      "agro_timecourse", "agro_multioutcome", "agro_inconsistent")) {
  name <- match.arg(name)
  map <- c(
    agro_nitrogen = "nitrogen", agro_irrigation = "irrigation", agro_biostimulants = "biostimulants",
    agro_crop_protection = "crop_protection", agro_soil_amendments = "soil_amendments",
    agro_multienv = "multienvironment", agro_timecourse = "timecourse",
    agro_multioutcome = "multioutcome", agro_inconsistent = "inconsistent"
  )
  nma_simulate(map[[name]], studies = if (name == "agro_multienv") 36 else 30,
               seed = 260819 + match(name, names(map)), inconsistency = if (name == "agro_inconsistent") 0.35 else 0)
}

#' Guided teaching tours
#'
#' @param topic Topic to tour.
#' @param run Run the lightweight data/preflight artifact.
#' @return A list of steps and optional artifact.
#' @export
nma_tour <- function(topic = c("network", "frequentist", "bayesian", "bootstrap", "inconsistency",
                               "dose", "component", "mlnmr", "publication"), run = FALSE) {
  topic <- match.arg(topic)
  steps <- switch(topic,
    network = c("Declare study/arm/treatment fields", "Audit duplicates and strata", "Inspect coarse/fine nodes",
                "Build geometry", "Diagnose disconnected components", "Apply only reviewed mappings", "Re-run preflight"),
    frequentist = c("Prepare effects", "Fit netmeta/NMA", "Inspect heterogeneity", "Assess inconsistency", "Estimate relative effects", "Report prediction intervals"),
    bayesian = c("Build multinma network", "Choose likelihood/link", "Justify priors", "Fit consistency model", "Check R-hat/ESS", "Assess model fit", "Node split where identifiable"),
    bootstrap = c("Choose study as resampling unit", "Set B and seed", "Refit complete studies", "Record failed/disconnected replicates", "Summarize coefficient/rank uncertainty"),
    inconsistency = c("Revisit transitivity", "Global design-by-treatment check", "Local node splitting", "Net heat", "Sensitivity to node policy"),
    dose = c("Preserve agent and dose", "Inspect dose overlap", "Fit simple shape", "Compare plausible shapes", "Check identifiability", "Plot uncertainty", "Restrict conclusions to dose support"),
    component = c("Define components before analysis", "Check identifiability", "Fit additive CNMA", "Add interactions only when supported", "Compare complex intervention predictions"),
    mlnmr = c("Harmonize IPD and AgD estimands", "Create IPD/AgD sources", "Combine network", "Specify effect modifiers", "Add integration", "Fit ML-NMR", "Check integration error and population-average effects"),
    publication = c("Network graph", "Forest/league table", "Transitivity/inconsistency figures", "Ranking uncertainty", "Sensitivity panel", "Reproducibility report")
  )
  artifact <- NULL
  if (isTRUE(run)) {
    dat <- nma_data_example("agro_irrigation")
    obj <- nma_data(dat, study = "study", treatment = "treatment", arm = "arm")
    artifact <- nma_preflight(obj, n = "n", mean = "mean", sd = "sd")
  }
  list(topic = topic, steps = steps, artifact = artifact)
}

#' Explain a methodological choice in concise teaching language
#' @param concept Concept keyword.
#' @return Character explanation.
#' @export
nma_explain <- function(concept = c("transitivity", "consistency", "multiarm", "ranking", "coarse_nodes", "bootstrap")) {
  concept <- match.arg(concept)
  switch(concept,
    transitivity = "Transitivity asks whether the studies supporting different treatment comparisons are sufficiently comparable in effect modifiers for indirect comparisons to be scientifically credible.",
    consistency = "Consistency is the statistical agreement between direct and indirect evidence after the transitivity assumption and model are accepted.",
    multiarm = "Arms from the same multi-arm study share randomization and controls; their treatment contrasts are correlated and must not be treated as independent studies.",
    ranking = "Ranking is a secondary summary. A high P-score or SUCRA does not by itself imply a large, precise, clinically or agronomically relevant treatment effect.",
    coarse_nodes = "Coarse nodes reduce unnecessary fragmentation for the primary synthesis, while fine nodes preserve dose, formulation, route or component detail for sensitivity and specialized models.",
    bootstrap = "For NMA, study-level bootstrap resamples complete studies so that all arms and correlated contrasts from a multi-arm trial remain together."
  )
}

#' Construct a recommended workflow object
#' @param x `nmaflow_data`.
#' @param question Optional analysis question.
#' @return Ordered workflow checklist with current preflight status.
#' @export
nma_workflow <- function(x, question = NULL) {
  x <- .nma_unwrap(x)
  pre <- nma_preflight(x)
  list(
    question = question,
    status = pre$status,
    steps = c("data architecture", "network formation", "effect measure", "transitivity", "primary model",
              "heterogeneity", "inconsistency", "sensitivity/bootstrap", "ranking with uncertainty",
              "publication outputs", "reproducibility"),
    preflight = pre
  )
}
