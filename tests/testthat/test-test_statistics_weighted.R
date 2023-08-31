test_that("Test statistics values (Informative sample) Model 1", {
  fit <- fit_facmod_pml(1, samp = "wtd", seed = 123)
  res <- all_tests(fit$fit, fit$svy)
  expect_snapshot(res)
})

# test_that("Test statistics values (Informative sample) Model 2", {
#   fit <- fit_facmod_pml(2, samp = "strcl", seed = 123)
#   res <- all_tests(fit$fit, fit$svy)
#   expect_snapshot(res)
# })

# test_that("Test statistics values (Informative sample) Model 3", {
#   fit <- fit_facmod_pml(3, samp = "strcl", seed = 123)
#   res <- all_tests(fit$fit, fit$svy)
#   expect_snapshot(res)
# })

# test_that("Test statistics values (Informative sample) Model 4", {
#   fit <- fit_facmod_pml(4, samp = "strcl", seed = 123)
#   res <- all_tests(fit$fit, fit$svy)
#   expect_snapshot(res)
# })

# test_that("Test statistics values (Informative sample) Model 5", {
#   fit <- fit_facmod_pml(5, samp = "strcl", seed = 123)
#   res <- all_tests(fit$fit, fit$svy)
#   expect_snapshot(res)
# })
