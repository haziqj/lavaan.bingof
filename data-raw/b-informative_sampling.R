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
library(kableExtra)

make_wt_res <- function(samp_size = 1000, type = "bias_se") {
  model_no <- 1
  dat <- gen_data_bin_wt(model_no, n = samp_size)
  fit0 <- lavaan::sem(model = txt_mod(model_no), data = dat, estimator = "PML",
                      std.lv = TRUE, sampling.weights = NULL)
  fit1 <- lavaan::sem(model = txt_mod(model_no), data = dat, estimator = "PML",
                      std.lv = TRUE, sampling.weights = "wt")
  ntilde <- lavaan.bingof:::neff(dat$wt)

  if (type == "bias_se") {
    # Bias and se results ------------------------------------------------------
    true_vals <- get_true_values(model_no, arrange = c("lambda", "tau", "rho"))
    se0 <- lavaan::parTable(fit0) %>% filter(free > 0)  # PML
    se1 <- lavaan::parTable(fit1) %>% filter(free > 0)  # PML (weighted)

    tmp <- tibble(
      name = names(lavaan::coef(fit0)),
      kind = names(true_vals) %>% gsub("\\d", "", .),
      truth = true_vals
    )
    res <- bind_rows(
      bind_cols(tmp, tibble(method = "PML",
                            est  = se0$est,
                            se  = se0$se)),
      bind_cols(tmp, tibble(method = "PMLW",
                            est  = se1$est,
                            se  = se1$se * sqrt(samp_size / ntilde)))
    ) %>%
      mutate(
        unc = qnorm(0.975) * se,
        covered = (truth >= est - se) & (truth <= est + se)
      )
  }
  if (type == "test_stats") {
    # Test statistics ----------------------------------------------------------
    res <- bind_rows(
      bind_cols(method = "PML", Pearson_test(fit0)),
      bind_cols(method = "PMLW", Pearson_test(fit1))
    )
  }

  as_tibble(res)

}
possfn <- possibly(make_wt_res, NA)

# Bias and se sims -------------------------------------------------------------
plan(multisession, workers = 30)
nsims <- 100

res <-
  future_map(seq_len(nsims), \(x) {
    possfn(samp_size = 1000, type = "bias_se")
  }, .progress = TRUE, .options = furrr_options(seed = NULL))

res_df <-
  res[sapply(res, is_tibble)] %>%
  do.call(rbind, .)

# Bias table
res_df %>%
  mutate(bias = est - truth) %>%
  group_by(method, name, kind, truth) %>%
  summarise(bias = mean((bias)), .groups = "drop") %>%
  pivot_wider(id_cols = c(name, kind, truth), names_from = method,
              values_from = bias)
  # kbl(format = "rst", digits = 2)

# SE table
res_df %>%
  group_by(method, name, kind, truth) %>%
  summarise(cov = mean(covered),
            sdse = sd(est) / mean(se),
            sd = sd(est),
            se = mean(se),
            .groups = "drop") %>%
  pivot_wider(id_cols = c(name, kind, truth), names_from = method,
              values_from = c(cov, sdse))

# Bias and se sims -------------------------------------------------------------
plan(multisession, workers = 30)
nsims <- 10

res <-
  future_map(seq_len(nsims), \(x) {
    possfn(samp_size = 1000, type = "test_stats")
  }, .progress = TRUE, .options = furrr_options(seed = NULL))

res_df <-
  res[sapply(res, is_tibble)] %>%
  do.call(rbind, .)

res <- list()
for (n in c(500, 1000)) {  #, 2500, 5000, 10000
  cat(paste("\nRunning with sample size n =", n, "\n"))

  out <-
    future_map(seq_len(nsims), \(x) {
      possfn(samp_size = n, type = "test_stats")
    }, .progress = TRUE, .options = furrr_options(seed = NULL))

  out <-
    out[sapply(out, is_tibble)] %>%
    do.call(rbind, .) %>%
    mutate(n = n)

  res <- c(res, list(out))
  cat("\n")
}


