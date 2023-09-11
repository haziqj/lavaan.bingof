# Simulations to investigate performance of weighted pairwise likelihood. We
# generate data according to an informative sampling scheme for Model 1 (5
# items). From this data, we fitted lavaan using method = "PML" (unweighted),
# and compare against the weighted version. To compare
#
# 1. Bias of parameters (true val vs PML bias vs PML weighted bias)
#
# 2. SE of paramaters (coverage and SD/SE) -- ratio of standard deviation of the
# parameter estimates in the simulation to the average standard errors.
#
# 3. Chi square rejection rates vs various sample sizes (n = 500, 1000, 2500,
# 5000, 10000)
#
# 4. Test statistics average value and df
library(tidyverse)
library(lavaan.bingof)
library(furrr)

make_test_stats_wt <- function(samp_size = 1000) {
  model_no <- 1
  dat <- gen_data_bin_wt(model_no, n = samp_size)
  fit0 <- lavaan::sem(model = txt_mod(model_no), data = dat, estimator = "PML",
                      std.lv = TRUE, sampling.weights = NULL)
  fit1 <- lavaan::sem(model = txt_mod(model_no), data = dat, estimator = "PML",
                      std.lv = TRUE, sampling.weights = "wt")

  # Bias results ---------------------------------------------------------------
  true_vals <- get_true_values(model_no, arrange = c("lambda", "tau", "rho"))
  coef0 <- lavaan::coef(fit0)  # PML
  coef1 <- lavaan::coef(fit1)  # PML (weighted)

  bias_df <- tibble(
    name = names(lavaan::coef(fit0)),
    kind = names(true_vals) %>% gsub("\\d", "", .),
    pl  = as.numeric(lavaan::coef(fit0)),
    plw = as.numeric(lavaan::coef(fit1)),
    truth = true_vals
  )

  bias_df
}; make_test_stats_wt()
