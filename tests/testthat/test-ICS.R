context("Information content of steps")

test_that("Information content of steps calculated correctly", {
  expect_equal(as.double(ICS(2, 2, 10000) * NUnrooted(4)), c(1, 2))
  expect_equal(signif(as.double(ICS(2, 3, 10000) * NUnrooted(5), 5)), c(2, 13))

})