test_that("reviewed node policy creates fine and coarse nodes", {
  d <- nma_simulate("nitrogen", studies = 8, seed = 11)
  x <- nma_data(d, study = "study", treatment = "treatment", arm = "arm")
  labs <- unique(d$treatment)
  dict <- data.frame(label = labs, fine = paste0(labs, "|fine"), coarse = labs)
  y <- nma_node_policy(x, dict, reviewed = TRUE, reason = "test dictionary")
  expect_true(all(c(".nma_node_fine", ".nma_node_coarse") %in% names(y$data)))
  expect_s3_class(nma_network(y, node = "coarse"), "nmaflow_network")
})

test_that("unit conversion is explicit and audited", {
  d <- nma_simulate("irrigation", studies = 6, seed = 12)
  d$unit <- "kg ha-1"
  x <- nma_data(d, study = "study", treatment = "treatment", arm = "arm", unit = "unit")
  expect_error(nma_convert_units(x, "mean", "unit", "kg ha-1", "t ha-1", .001), "reviewed")
  y <- nma_convert_units(x, "mean", "unit", "kg ha-1", "t ha-1", .001,
                         reviewed = TRUE, reason = "test conversion")
  expect_true(all(y$data$unit == "t ha-1"))
  expect_true(any(y$audit$action == "unit_conversion"))
})

test_that("strata remain separate", {
  d <- nma_simulate("multioutcome", studies = 6, seed = 13)
  x <- nma_data(d, study = "study", treatment = "treatment", arm = "arm",
                outcome = "outcome", unit = "unit")
  ss <- nma_split_strata(x, by = c("outcome", "unit"))
  expect_gte(length(ss), 3)
  expect_true(all(vapply(ss, inherits, logical(1), "nmaflow_data")))
})

test_that("repair report compares geometry", {
  d <- nma_simulate("irrigation", studies = 6, seed = 14)
  x <- nma_data(d, study = "study", treatment = "treatment", arm = "arm")
  r <- nma_network_repair_report(x, x)
  expect_s3_class(r, "nmaflow_repair_report")
  expect_equal(r$before$n_nodes, r$after$n_nodes)
})
