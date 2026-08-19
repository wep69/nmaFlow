test_that("simulation is reproducible", {
  a <- nma_simulate("nitrogen", studies = 8, seed = 100)
  b <- nma_simulate("nitrogen", studies = 8, seed = 100)
  expect_equal(a, b)
})

test_that("teaching examples are available", {
  x <- nma_data_example("agro_biostimulants")
  expect_true(is.data.frame(x))
  expect_true(all(c("study", "treatment") %in% names(x)))
})
