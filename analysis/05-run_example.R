# This script to do single runs (before starting the simulations)
source("01-utilities.R")
source("02-data_gen_srs.R")
source("03-data_gen_complex.R")
source("04-test_statistics.R")

run_single_fit <- function(i, seed = 4423, H1 = FALSE) {
  mod <- txt_mod(i)

  # Simple random sampling
  dat <- gen_data_bin(i, seed = seed, H1 = H1)
  fit <- sem(model = mod, data = dat, estimator = "PML", std.lv = TRUE,
             sampling.weights = NULL)
  res1 <- all_tests(fit)

  # Create population for complex sampling
  pop <- make_population(i, H1 = H1)

  # Stratified sampling
  dat <- gen_data_bin_complex1(pop, seed = seed)
  fit <- sem(model = mod, data = dat, estimator = "PML", std.lv = TRUE,
             sampling.weights = "wt")
  res2 <- all_tests(fit)

  # Cluster sampling
  dat <- gen_data_bin_complex2(pop, seed = seed)
  fit <- sem(model = mod, data = dat, estimator = "PML", std.lv = TRUE,
             sampling.weights = "wt")
  res3 <- all_tests(fit)

  # Stratified-cluster sampling
  dat <- gen_data_bin_complex3(pop, seed = seed)
  fit <- sem(model = mod, data = dat, estimator = "PML", std.lv = TRUE,
             sampling.weights = "wt")
  res4 <- all_tests(fit)

  list(srs = res1, strat  = res2, clust = res3, strclu = res4)
}

run_single_fit(1)
run_single_fit(1, H1 = TRUE, seed = 4423)

