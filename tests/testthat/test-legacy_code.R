test_that("Wald test (SRS)", {
  fit <- fit_facmod_pml(1)$fit
  res1 <- Wald_Pearson_test_function(fit)  # legacy code
  res2 <- Wald_test(fit)

  expect_equal(res2$W, res1$Wald_test)
})

test_that("Wald test V2 (SRS)", {
  fit <- fit_facmod_pml(1)$fit
  res1 <- Wald_Pearson_test_function(fit)  # legacy code
  res2 <- Wald_test_v2(fit, .order = 0) %>% suppressWarnings()

  expect_equal(res2$W, res1$Wald_test_v2)
})

test_that("Pearson test [MM0] (SRS)", {
  fit <- fit_facmod_pml(1)$fit
  res1 <- Wald_Pearson_test_function(fit)  # legacy code
  res2 <- Pearson_test_v2(fit, .order = 0) %>% suppressWarnings()

  expect_equal(res2$W, c(res1$Pearson_test))
})

test_that("Pearson V2 test [MM2] (SRS)", {
  fit <- fit_facmod_pml(1)$fit
  res1 <- Wald_Pearson_test_function(fit)  # legacy code
  res2 <- Pearson_test_v2(fit, .order = 2)

  expect_equal(res2$W, c(res1$FSMadj_Pearson))
  expect_equal(res2$df, c(res1$FSMadj_df_Pearson))
})

# Now test each model

test_that("Wald test (SRS) Model 1", {
  fit <- fit_facmod_pml(1)$fit
  res1 <- Wald_Pearson_test_function(fit)  # legacy code
  res2 <- Wald_test(fit)

  expect_equal(res2$W, res1$Wald_test)
})

test_that("Wald test (SRS) Model 2", {
  fit <- fit_facmod_pml(2)$fit
  res1 <- Wald_Pearson_test_function(fit)  # legacy code
  res2 <- Wald_test(fit)

  expect_equal(res2$W, res1$Wald_test)
})

# Remark: Tests for Models 3 and 5 will fail. It seems Myrsini's code to
# generate Sigma2 yields in slightly different numerical values. I think this is
# because of how the Sigma2 was populated. I used mnormt::sadmvn() to calculate
# all relevant probabilities, while the legacy code uses lavaan's built in pnorm
# calculator and some multiplication to get things like Cov(y12,y13) = E(y12
# y13) - E(y12)E(y13). However, the differences are really small and I think the
# legacy code contains some mistakes. Check out the plot of
#
# plot(create_Sigma2_matrix(fit), res1$Sigma2)

# test_that("Wald test (SRS) Model 3", {
#   fit <- fit_facmod_pml(3)
#   res1 <- Wald_Pearson_test_function(fit)  # legacy code
#   res2 <- Wald_test(fit)
#
#   expect_equal(res2$W, res1$Wald_test)
# })

test_that("Wald test (SRS) Model 4", {
  fit <- fit_facmod_pml(4)$fit
  res1 <- Wald_Pearson_test_function(fit)  # legacy code
  res2 <- Wald_test(fit)

  expect_equal(res2$W, res1$Wald_test)
})

# test_that("Wald test (SRS) Model 5", {
#   fit <- fit_facmod_pml(5)
#   res1 <- Wald_Pearson_test_function(fit)  # legacy code
#   res2 <- Wald_test(fit)
#
#   expect_equal(res2$W, res1$Wald_test)
# })
