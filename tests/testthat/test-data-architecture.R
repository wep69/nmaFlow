test_that("basic agronomic network passes structural creation", {
  d <- nma_simulate("irrigation", studies = 12, seed = 1)
  x <- nma_data(d, study = "study", treatment = "treatment", arm = "arm")
  expect_s3_class(x, "nmaflow_data")
  expect_true(nrow(x$data) >= 24)
  expect_s3_class(nma_network(x), "nmaflow_network")
})

test_that("exact duplicate rows are quarantined", {
  d <- nma_simulate("irrigation", studies = 8, seed = 2)
  d <- rbind(d, d[1, ])
  x <- nma_data(d, study = "study", treatment = "treatment", arm = "arm")
  y <- nma_repair(x)
  expect_equal(nrow(y$quarantine), 1)
  expect_equal(nrow(y$data), nrow(d) - 1)
})

test_that("semantic mapping requires reviewed guard", {
  d <- nma_simulate("irrigation", studies = 8, seed = 3)
  x <- nma_data(d, study = "study", treatment = "treatment", arm = "arm")
  expect_error(nma_apply_mapping(x, c("Full irrigation" = "Irrigated")), "reviewed")
  y <- nma_apply_mapping(x, c("Full irrigation" = "Irrigated"), reviewed = TRUE, reason = "protocol node policy")
  expect_true(any(y$data$treatment == "Irrigated"))
})

test_that("disconnected networks are not auto-reconnected", {
  d <- data.frame(study = c("S1","S1","S2","S2"), arm = c(1,2,1,2),
                  treatment = c("A","B","C","D"), mean = c(1,2,1,2), sd = 1, n = 20)
  x <- nma_data(d, study = "study", treatment = "treatment", arm = "arm")
  net <- nma_network(x)
  expect_false(net$connected)
  plan <- nma_reconnect_plan(x)
  expect_false(plan$connected)
  expect_equal(nrow(x$data), 4)
})
