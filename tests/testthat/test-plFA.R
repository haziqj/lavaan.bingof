suppressPackageStartupMessages(library(plFA))
dat <- gen_data_bin(n = 1000, seed = 123)
mod <- "eta =~ y1 + y2 + y3 + y4 + y5"
fit1 <- plFA::cfa(mod, dat, std.lv = TRUE)
fit2 <- lavaan::cfa(mod, dat, std.lv = TRUE, estimator = "PML")

test_that("Coefficients and implied cov mat similar", {
  expect_equal(lavaan::coef(fit1), lavaan::coef(fit2), tolerance = 1e-5)
  expect_equal(fit1@Fit@Sigma.hat, fit2@Fit@Sigma.hat, tolerance = 1e-5)
})

test_that("Extracted lavaan info similar", {
  x <- extract_lavaan_info(fit1)
  y <- extract_lavaan_info(fit2)

  # No need to check the following
  x$lavoptions <- x$lavpartable <- NULL
  y$lavoptions <- y$lavpartable <- NULL

  expect_equal(x, y, tolerance = 1e-5)
})

test_that("All tests give similar output", {
  these_tests <- c(
    "Wald_test",
    "Wald_vcovf_test",
    "Wald_diag_test",
    "Pearson_test",
    "RSS_test",
    "Multn_test"
  )
  names(these_tests) <- these_tests
  for (i in seq_along(these_tests)) {
    x <- these_tests[i]
    out1 <- do.call(x, list(fit1))
    out2 <- do.call(x, list(fit2))
    expect_equal(out1, out2, tolerance = 1e-3)
  }
})

test_that("Similar calc test stuff output", {
  x <- calc_test_stuff(fit1)
  y <- calc_test_stuff(fit2)
  expect_equal(x, y, tolerance = 1e-4)
})
