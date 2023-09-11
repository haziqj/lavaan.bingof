test_that("Test statistics values (Stratified sampling v2) Model 1", {
  fit <- fit_facmod_pml(1, n = 500, samp = "strat2", seed = 123)
  res <- all_tests(fit$fit)
  expect_snapshot(res)
})
