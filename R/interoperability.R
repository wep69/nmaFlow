#' Export data and analysis templates to external ecosystems
#'
#' @param x `nmaflow_data`, `nmaflow_effects`, or data frame.
#' @param engine `python`, `julia`, `sas`, `stata`, `mplus`, `jags`, or `bugs`.
#' @param dir Output directory.
#' @param prefix File prefix.
#' @return Paths created. Templates are starting points and must be scientifically reviewed.
#' @export
nma_export <- function(x, engine = c("python", "julia", "sas", "stata", "mplus", "jags", "bugs"),
                       dir = "nma_export", prefix = "nmaflow") {
  engine <- match.arg(engine); dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  d <- if (inherits(x, "nmaflow_data")) x$data else if (inherits(x, "nmaflow_effects")) x$data else x
  .nma_assert_data(d)
  csv <- file.path(dir, paste0(prefix, "_data.csv")); utils::write.csv(d, csv, row.names = FALSE, na = "")
  script <- switch(engine,
    python = c(
      "# nmaFlow Python interoperability template", "# Python is optional; review likelihood, priors and multi-arm structure before fitting.",
      "import pandas as pd", "import pymc as pm", sprintf("dat = pd.read_csv(r'%s')", normalizePath(csv, winslash = "/", mustWork = FALSE)),
      "# Build the contrast/design matrices explicitly before defining a hierarchical PyMC model.",
      "# Preserve shared-control covariance for multi-arm studies."
    ),
    julia = c(
      "# nmaFlow Julia interoperability template", "# Julia is optional; this template is intended for Turing/JuliaBUGS model development.",
      "using CSV, DataFrames", sprintf("dat = CSV.read(raw\"%s\", DataFrame)", normalizePath(csv, winslash = "/", mustWork = FALSE)),
      "# Define treatment coding, multi-arm covariance, likelihood and priors before @model specification."
    ),
    sas = c(
      "/* nmaFlow SAS template */", "/* Review coding and covariance before PROC BGLIMM. */",
      sprintf("proc import datafile=\"%s\" out=nma dbms=csv replace; guessingrows=max; run;", normalizePath(csv, winslash = "/", mustWork = FALSE)),
      "/* Construct study-specific treatment contrasts and a Bayesian hierarchical model using PROC BGLIMM. */"
    ),
    stata = c(
      "* nmaFlow Stata interoperability template", sprintf("import delimited using \"%s\", clear", normalizePath(csv, winslash = "/", mustWork = FALSE)),
      "* Declare the network using the installed/official Stata meta-analysis workflow appropriate to your data format.",
      "* Verify treatment coding, effect direction, multi-arm handling, and inconsistency diagnostics."
    ),
    mplus = c(
      "TITLE: nmaFlow SEM-NMA interoperability template;", sprintf("DATA: FILE = %s;", basename(csv)),
      "VARIABLE: ! Define study-level contrasts and covariance variables here;",
      "ANALYSIS: ESTIMATOR = MLR;", "MODEL: ! Specify the SEM-NMA path structure after validating coding;",
      "OUTPUT: TECH1 TECH8 CINTERVAL;"
    ),
    jags = c(
      "# nmaFlow JAGS template", "# This file does not invent a likelihood. Fill in data-specific NMA equations after design review.",
      "model {", "  # study baselines", "  # basic treatment effects", "  # multi-arm correction", "  # likelihood and link", "}"
    ),
    bugs = c(
      "# nmaFlow OpenBUGS/WinBUGS template", "model {", "  # Define arm-level likelihood and consistency equations.",
      "  # Add multi-arm corrections and priors explicitly.", "}"
    )
  )
  ext <- c(python = "py", julia = "jl", sas = "sas", stata = "do", mplus = "inp", jags = "jags", bugs = "bug")[[engine]]
  path <- file.path(dir, paste0(prefix, "_", engine, ".", ext)); writeLines(script, path)
  c(data = normalizePath(csv, mustWork = FALSE), script = normalizePath(path, mustWork = FALSE))
}

#' Save a reproducibility bundle
#'
#' @param x nmaFlow object.
#' @param dir Bundle directory.
#' @param include_raw Include raw imported data when present.
#' @return Paths created.
#' @export
nma_reproduce <- function(x, dir = "nmaflow_repro", include_raw = TRUE) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  obj <- file.path(dir, "analysis_object.rds"); saveRDS(x, obj)
  paths <- c(object = obj)
  if (inherits(x, "nmaflow_data")) {
    utils::write.csv(x$data, file.path(dir, "analysis_data.csv"), row.names = FALSE, na = "")
    utils::write.csv(x$audit, file.path(dir, "audit_log.csv"), row.names = FALSE, na = "")
    utils::write.csv(x$quarantine, file.path(dir, "quarantine.csv"), row.names = FALSE, na = "")
    if (isTRUE(include_raw)) utils::write.csv(x$raw_data, file.path(dir, "raw_import.csv"), row.names = FALSE, na = "")
  }
  writeLines(utils::capture.output(nma_session()), file.path(dir, "session.txt"))
  manifest <- c(
    paste0("created_utc: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
    paste0("class: ", paste(class(x), collapse = ", ")),
    paste0("nmaFlow_version: ", tryCatch(as.character(utils::packageVersion("nmaFlow")), error = function(e) "development-tree"))
  )
  writeLines(manifest, file.path(dir, "manifest.txt"))
  normalizePath(dir, mustWork = FALSE)
}

#' Reproducibility session report
#' @return A list containing R session, nmaFlow options and optional backend versions.
#' @export
nma_session <- function() {
  caps <- nma_capabilities(); available <- caps$package[caps$available]
  versions <- stats::setNames(vapply(available, function(p) as.character(utils::packageVersion(p)), character(1)), available)
  list(session = utils::sessionInfo(), options = options()[c("nmaFlow.digits", "nmaFlow.seed", "nmaFlow.auto_repair")],
       backend_versions = versions)
}
