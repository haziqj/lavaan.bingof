test_that("Test statistics values (Stratified sampling) Model 1", {
  fit <- fit_facmod_pml(1, samp = "strat", seed = 123)
  res <- all_tests(fit$fit)
  expect_snapshot(res)
})

test_that("Test statistics values (Stratified sampling) Model 2", {
  fit <- fit_facmod_pml(2, samp = "strat", seed = 123)
  res <- all_tests(fit$fit)
  expect_snapshot(res)
})

# test_that("Test statistics values (Stratified sampling) Model 3", {
#   fit <- fit_facmod_pml(3, samp = "strat", seed = 123)
#   res <- all_tests(fit$fit)
#   expect_snapshot(res)
# })

test_that("Test statistics values (Stratified sampling) Model 4", {
  fit <- fit_facmod_pml(4, samp = "strat", seed = 123)
  res <- all_tests(fit$fit)
  expect_snapshot(res)
})

# test_that("Test statistics values (Stratified sampling) Model 5", {
#   fit <- fit_facmod_pml(5, samp = "strat", seed = 123)
#   res <- all_tests(fit$fit)
#   expect_snapshot(res)
# })
