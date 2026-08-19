test_that("80-block registry is complete and unique", {
  b <- nma_blocks()
  expect_equal(nrow(b), 80)
  expect_equal(length(unique(b$block)), 80)
  expect_false(anyDuplicated(b[["function"]]) > 0)
})

test_that("capability table includes key backends", {
  c <- nma_capabilities()
  expect_true(all(c("netmeta", "multinma", "MBNMAdose", "crossnma") %in% c$package))
})
