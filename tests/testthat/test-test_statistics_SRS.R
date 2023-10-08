# test_that("Test statistics values (SRS) Model 1", {
#   fit <- fit_facmod_pml(1, seed = 123)
#   res <- all_tests(fit$fit)
#   expect_snapshot(res)
# })
#
# test_that("Test statistics values (SRS) Model 2", {
#   fit <- fit_facmod_pml(2, seed = 123)
#   res <- all_tests(fit$fit)
#   expect_snapshot(res)
# })

test_that("Test statistics values (SRS) Model 3", {
  fit <- fit_facmod_pml(3, seed = 123)
  res <- all_tests(fit$fit)
  expect_snapshot(res)
})

# test_that("Test statistics values (SRS) Model 4", {
#   fit <- fit_facmod_pml(4, seed = 123)
#   res <- all_tests(fit$fit)
#   expect_snapshot(res)
# })

test_that("Test statistics values (SRS) Model 5", {
  fit <- fit_facmod_pml(5, seed = 123)
  res <- all_tests(fit$fit)
  expect_snapshot(res)
})
