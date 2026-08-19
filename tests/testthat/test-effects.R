test_that("continuous arm summaries produce MD contrasts", {
  d <- nma_simulate("irrigation", studies = 10, seed = 4)
  x <- nma_data(d, study = "study", treatment = "treatment", arm = "arm")
  e <- nma_effects(x, "MD", mean = "mean", sd = "sd", n = "n")
  expect_s3_class(e, "nmaflow_effects")
  expect_true(all(c("TE", "seTE", "treat1", "treat2") %in% names(e$data)))
  expect_true(all(e$data$seTE > 0, na.rm = TRUE))
})

test_that("dispersion recovery follows explicit hierarchy", {
  expect_equal(nma_dispersion(25, se = 0.2)$sd, 1)
  expect_equal(nma_dispersion(25, mse = 4)$sd, 2)
  expect_equal(nma_dispersion(25)$method, "unavailable")
})
