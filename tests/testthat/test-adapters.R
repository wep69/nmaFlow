test_that("nma_fit_nma maps arm-level columns to NMA setup", {
  skip_if_not_installed("NMA")
  d <- nma_data_example("agro_irrigation")
  x <- nma_data(d, study = "study", treatment = "treatment", arm = "arm")
  fit <- nma_fit_nma(
    x, measure = "MD", ref = "Rainfed",
    mean = "mean", sd = "sd", n = "n"
  )
  expect_s3_class(fit, "nmaflow_fit")
  expect_identical(fit$engine, "NMA")
})

test_that("nma_dose maps contrast columns to netdose", {
  skip_if_not_installed("netdose")
  doses <- c(0, 25, 50, 100, 200)
  d <- expand.grid(
    study = sprintf("S%02d", 1:8),
    treatment = factor(doses, levels = doses),
    stringsAsFactors = FALSE
  )
  d$arm <- rep(seq_along(doses), times = 8)
  d$mean <- 2 + 0.8 * as.numeric(as.character(d$treatment)) / 200
  d$sd <- 0.25
  d$n <- 30
  obj <- nma_data(d, study = "study", treatment = "treatment", arm = "arm")
  eff <- nma_effects(obj, measure = "MD", mean = "mean", sd = "sd", n = "n")
  con <- eff$data
  con$dose1 <- as.numeric(as.character(con$treat1))
  con$dose2 <- as.numeric(as.character(con$treat2))
  con$agent1 <- ifelse(con$dose1 == 0, "Control", "Dose")
  con$agent2 <- ifelse(con$dose2 == 0, "Control", "Dose")
  con$studlab <- con$study
  con$sm <- "MD"
  fit <- nma_dose(
    con[, c("studlab", "agent1", "dose1", "agent2", "dose2", "TE", "seTE", "sm")],
    framework = "frequentist", model = "linear", common = TRUE, random = FALSE
  )
  expect_s3_class(fit, "nmaflow_fit")
  expect_identical(fit$engine, "netdose")
})

test_that("nma_time does not pass unsupported method", {
  skip_if_not_installed("MBNMAtime")
  d <- expand.grid(
    studyID = c("A", "B"),
    treatment = c("T1", "T2"),
    time = c(0, 7, 14, 21, 28),
    stringsAsFactors = FALSE
  )
  d$y <- 2 + 0.04 * d$time + (d$treatment == "T2") * 0.1
  d$se <- 0.2
  d$n <- 20
  fit <- nma_time(
    d, fun = MBNMAtime::tpoly(degree = 1),
    n.iter = 40, n.burnin = 20, n.chains = 2
  )
  expect_s3_class(fit, "nmaflow_fit")
  expect_identical(fit$engine, "MBNMAtime")
})

test_that("nma_crossdesign preserves named backend arguments", {
  skip_if_not_installed("crossnma")
  ipd_env <- new.env(parent = emptyenv())
  std_env <- new.env(parent = emptyenv())
  utils::data("ipddata", package = "crossnma", envir = ipd_env)
  utils::data("stddata", package = "crossnma", envir = std_env)
  ipddata <- ipd_env$ipddata
  stddata <- std_env$stddata
  model <- nma_crossdesign(
    trt = "treat", study = "id", outcome = "relapse", n = "n", design = "design",
    prt.data = ipddata, std.data = stddata,
    method.bias = "naive", sm = "OR", run = FALSE
  )
  expect_true(inherits(model, "crossnma.model") || is.list(model))
})

test_that("nma_table extracts netmeta ranking scores", {
  d <- nma_data_example("agro_irrigation")
  x <- nma_data(d, study = "study", treatment = "treatment", arm = "arm")
  e <- nma_effects(x, measure = "MD", mean = "mean", sd = "sd", n = "n")
  fit <- nma_fit_freq(e, sm = "MD", random = TRUE)
  tab <- nma_table(fit, component = "rank", format = "data.frame")
  expect_s3_class(tab, "data.frame")
  expect_true("P.score" %in% names(tab))
})

test_that("nma_threshold identifies nma_thresh and validates inputs", {
  skip_if_not_installed("nmathresh")
  expect_error(
    nma_threshold(structure(list(), class = "nmaflow_fit")),
    "mean.dk"
  )
})

test_that("nma_predict reports unsupported netmeta backend clearly", {
  d <- nma_data_example("agro_irrigation")
  x <- nma_data(d, study = "study", treatment = "treatment", arm = "arm")
  e <- nma_effects(x, measure = "MD", mean = "mean", sd = "sd", n = "n")
  fit <- nma_fit_freq(e, sm = "MD", random = TRUE)
  expect_error(nma_predict(fit), "Prediction is not available for netmeta")
})

test_that("nma_multivariate supports vector y through a formula adapter", {
  skip_if_not_installed("mixmeta")
  y <- c(0.1, 0.2, 0.15, 0.08)
  fit <- nma_multivariate(y = y, S = rep(0.01, length(y)))
  expect_s3_class(fit, "nmaflow_fit")
  expect_identical(fit$engine, "mixmeta")
})

test_that("backend result checker fails loudly", {
  expect_error(
    nmaFlow:::.nma_check_backend_result(list(engine = "failed to create sampler"), "test backend"),
    "failed to create sampler"
  )
})
