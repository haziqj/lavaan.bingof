test_that("Frequency count", {
  dat <- gen_data_bin(1)
  tab1 <- create_pairwise_table(dat)  # bingof
  tab2 <- lavaan::lavTables(dat)  # lavaan
  expect_equal(tab1$freq, tab2$obs.freq)
})

# test_that("Weighted frequency count", {
#   dat <- gen_data_bin_wt(1)
#   tab1 <- create_pairwise_table(select(dat, -wt), dat$wt)  # bingof
#   tab2 <- lavaan::lavTables(dat)  # lavaan
#   expect_equal(tab1$freq, tab2$obs.freq)
# })

test_that("Sensitivity matrix is the same", {
  model_no <- 1
  samp_size <- 1000
  mod <- fit_facmod_pml(model_no, n = samp_size)

  Hinv_lav <- get_sensitivity_inv_mat(mod$fit)
  Hinv_pml <- get_Hinv_mat(select(mod$dat, starts_with("y")), model_no, NULL,
                           .trace = FALSE)


  expect_equal(diag(Hinv_lav), diag(Hinv_pml, names = FALSE), tolerance = 1e-3)
})
