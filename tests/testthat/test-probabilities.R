# Check sample/population prevalences match the theoretical probabilities

test_that("Proportions align (SRS)", {
  dat <- gen_data_bin(5, n = 200000)
  tau <- get_true_values(5)
  tau <- tau[grepl("tau", names(tau))]

  true_prop <- pnorm(tau, lower.tail = FALSE)
  obs_prop  <- apply(dat == 1, 2, mean)
  names(true_prop) <- names(obs_prop)
  expect_equal(obs_prop, true_prop, tolerance = 0.01)
})

test_that("Tetrachoric correlations align (SRS)", {
  dat <- gen_data_bin(5, n = 500000, return_all = TRUE)
  mat <- dat %>% select(starts_with("ystar")) %>% cor()
  Sigmay <- get_Sigmay(5)
  expect_equal(mat[upper.tri(mat)], Sigmay[upper.tri(Sigmay)],
               tolerance = 0.05)
})

test_that("Proportions align (Population)", {
  dat <- make_population(1)
  tau <- get_tau(1)

  # proportion zeroes
  true_prop <- pnorm(tau)
  obs_prop  <- dat %>%
    select(starts_with("y")) %>%
    mutate(across(everything(), \(x) x == 0)) %>%
    apply(2, mean)
  names(true_prop) <- names(obs_prop)
  expect_equal(obs_prop, true_prop, tolerance = 0.01)
})
