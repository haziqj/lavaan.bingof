library(tidyverse)
library(lavaan.bingof)
library(furrr)
library(cli)
library(gt)
plan(multisession, workers = parallel::detectCores() - 1)
nsims <- 1000

# 1. Run additional simulations for smaller sample sizes n = 500 and n =
# 1000(with the n = 5000). Table showing mean bias, mean abs. bias, mean
# coverage, mean sd/se ratio.

# 2. Another table showing Wald and Pearson (unadjusted vs weighted)

make_wt_res <- function(samp_size = 5000, type = "bias_se", .popfac = 1 / 0.01,
                        model_no = 1) {
  dat <- gen_data_bin_wt(model_no, n = samp_size, popfac = .popfac)
  fit0 <- lavaan::sem(model = txt_mod(model_no), data = dat, estimator = "PML",
                      std.lv = TRUE, sampling.weights = NULL)
  fit1 <- lavaan::sem(model = txt_mod(model_no), data = dat, estimator = "PML",
                      std.lv = TRUE, sampling.weights = "wt")

  if (type == "bias_se") {
    # Bias and se results ------------------------------------------------------
    true_vals <- get_true_values(model_no, arrange = c("lambda", "tau", "rho"))
    se0 <- lavaan::parTable(fit0) |> filter(free > 0)  # PML
    se1 <- lavaan::parTable(fit1) |> filter(free > 0)  # PML (weighted)

    tab <- tibble(
      name  = names(lavaan::coef(fit0)),
      kind  = gsub("\\d", "", names(true_vals)),
      truth = true_vals,
      n     = samp_size
    )
    res <-
      bind_rows(
        bind_cols(tab, tibble(method = "PML",
                              est    = se0$est,
                              se     = se0$se)),
        bind_cols(tab, tibble(method = "PMLW",
                              est    = se1$est,
                              se     = se1$se))
      ) |>
      mutate(
        unc = qnorm(0.975) * se,
        covered = (truth >= est - unc) & (truth <= est + unc)
      )
  }
  if (type == "test_stats") {
    # Test statistics ----------------------------------------------------------
    # PML
    tab0 <- all_tests(fit0)
    tab1 <- all_tests(fit1, Sigma2 = "weighted")

    res <- bind_rows(
      bind_cols(n = samp_size, method = "PML", tab0),
      bind_cols(n = samp_size, method = "PMLW", tab1),
    )
  }

  # Check for Heywood cases
  res$var_ok <-
    lavaan:::lav_object_post_check(fit0) & lavaan:::lav_object_post_check(fit1)
  res$converged <-
    fit0@Fit@converged & fit1@Fit@converged
  as_tibble(res)

}
make_wt_res <- possibly(make_wt_res, NA)

# Run simulations --------------------------------------------------------------

# Bias, coverage and SE
res <- list()
for (n in c(500, 1000, 5000)) {
  cli_alert_info("Running with sample size n = {n}")

  out <-
    future_map(seq_len(nsims), \(x) {
      make_wt_res(samp_size = n, type = "bias_se", .popfac = 500000 / n,
                  model_no = 5)
    }, .progress = TRUE, .options = furrr_options(seed = NULL))

  out <-
    out[sapply(out, is_tibble)] |>
    do.call(what = rbind)

  res <- c(res, list(out))
  cat("\n")
}
res_bias_se_df <- do.call(rbind, res)

# Test statistics
res <- list()
for (n in c(500, 1000, 5000)) {
  cli_alert_info("Running with sample size n = {n}")

  out <-
    future_map(seq_len(nsims), \(x) {
      make_wt_res(samp_size = n, type = "test_stats", .popfac = 500000 / n)
    }, .progress = TRUE, .options = furrr_options(seed = NULL))

  out <-
    out[sapply(out, is_tibble)] |>
    do.call(what = rbind)

  res <- c(res, list(out))
  cat("\n")
}
res_test_stats <- do.call(rbind, res)

# Create bias table ------------------------------------------------------------
tab_bias <-
  res_bias_se_df |>
  summarise(
    truth = first(truth),
    bias = mean(est - truth),
    # abs_bias = mean(abs(est - truth)),
    .by = c(n, method, name)
  ) |>
  pivot_wider(names_from = c(method, n), values_from = c(bias))

tab_bias_gt <-
  tab_bias |>
  filter(
    name %in% c("eta1=~y1", "eta1=~y2", "eta1=~y3", "eta1=~y4", "eta1=~y5",
                "y1|t1", "y2|t1", "y3|t1", "y4|t1", "y5|t1",
                "eta1~~eta2", "eta1~~eta3", "eta2~~eta3")
  ) |>
  mutate(
    name = gsub("eta1=~y", "$\\\\lambda_", name),
    name = gsub("y", "$\\\\tau_", name),
    name = gsub("\\|t1", "", name),
    name = gsub("eta1~~eta2$", "$\\\\rho_{12}", name),
    name = gsub("eta1~~eta3$", "$\\\\rho_{13}", name),
    name = gsub("eta2~~eta3$", "$\\\\rho_{23}", name),
    name = paste0(name, "$")
  ) |>
  gt(rowname_col = "name") |>
  fmt_number(columns = -"truth", decimals = 3) |>
  fmt_markdown("name") |>
  tab_spanner(
    label = md("$n = 500$"),
    columns = ends_with("_500")
  ) |>
  tab_spanner(
    label = md("$n = 1000$"),
    columns = ends_with("_1000")
  ) |>
  tab_spanner(
    label = md("$n = 5000$"),
    columns = ends_with("_5000")
  ) |>
  tab_spanner(
    label = "Abs. Bias",
    columns = starts_with("abs_bias")
  ) |>
  tab_spanner(
    label = "Bias",
    columns = starts_with("bias")
  ) |>
  tab_row_group(
    label = md("**Factor correlations**"),
    rows = contains("rho")
  ) |>
  tab_row_group(
    label = md("**Thresholds**"),
    rows = contains("tau")
  ) |>
  tab_row_group(
    label = md("**Loadings**"),
    rows = contains("lambda")
  ) |>
  # cols_width(truth ~ px(70)) |>
  cols_label(
    contains("PML_") ~ "PML",
    contains("PMLW_") ~ "PMLW",
    truth = "True values"
  )

# Create se table --------------------------------------------------------------
tab_se <-
  res_bias_se_df |>
  summarise(
    sd = sqrt(mean((est - truth) ^ 2)),
    se = mean(se),
    cov = mean(covered),
    .by = c(n, method, name)
  ) |>
  mutate(ratio = sd / se) |>
  select(-sd, -se) |>
  pivot_wider(names_from = c(method, n), values_from = c(cov, ratio))

tab_se_gt <-
  tab_se |>
  filter(
    name %in% c("eta1=~y1", "eta1=~y2", "eta1=~y3", "eta1=~y4", "eta1=~y5",
                "y1|t1", "y2|t1", "y3|t1", "y4|t1", "y5|t1",
                "eta1~~eta2", "eta1~~eta3", "eta2~~eta3")
  ) |>
  mutate(
    name = gsub("eta1=~y", "$\\\\lambda_", name),
    name = gsub("y", "$\\\\tau_", name),
    name = gsub("\\|t1", "", name),
    name = gsub("eta1~~eta2$", "$\\\\rho_{12}", name),
    name = gsub("eta1~~eta3$", "$\\\\rho_{13}", name),
    name = gsub("eta2~~eta3$", "$\\\\rho_{23}", name),
    name = paste0(name, "$")
  ) |>
  gt(rowname_col = "name") |>
  fmt_number(decimals = 2) |>
  fmt_markdown("name") |>
  tab_spanner(
    label = md("$n = 500$"),
    columns = ends_with("_500")
  ) |>
  tab_spanner(
    label = md("$n = 1000$"),
    columns = ends_with("_1000")
  ) |>
  tab_spanner(
    label = md("$n = 5000$"),
    columns = ends_with("_5000")
  ) |>
  tab_spanner(
    label = "Coverage",
    columns = starts_with("cov_")
  ) |>
  tab_spanner(
    label = "SD/SE",
    columns = starts_with("ratio_")
  ) |>
  tab_row_group(
    label = md("**Factor correlations**"),
    rows = contains("rho")
  ) |>
  tab_row_group(
    label = md("**Thresholds**"),
    rows = contains("tau")
  ) |>
  tab_row_group(
    label = md("**Loadings**"),
    rows = contains("lambda")
  ) |>
  cols_label(
    contains("PML_") ~ "PML",
    contains("PMLW_") ~ "PMLW"
  )

# Create test stats table ------------------------------------------------------
tab_test_stats <-
  res_test_stats |>
  filter(name %in% c("WaldVCF", "Pearson,MM3")) |>
  mutate(name = gsub(",MM3", "", name)) |>
  summarise(
    rejrate = mean(pval < 0.05),
    X2 = mean(X2),
    df = mean(df),
    .by = c(n, method, name)
  ) |>
  mutate(X2df = glue::glue("{iprior::dec_plac(X2,2)} ({iprior::dec_plac(df,2)})")) |>
  select(-X2, -df) |>
  pivot_wider(names_from = c(name, method), values_from = c(rejrate, X2df))

tab_test_stats_gt <-
  gt(tab_test_stats) |>
  # fmt_auto() |>
  fmt_number(columns = -n, decimals = 3) |>
  tab_spanner(
    label = "Wald",
    columns = contains("Wald")
  ) |>
  tab_spanner(
    label = "Pearson",
    columns = contains("Pearson")
  ) |>
  tab_spanner(
    label = "Rej. Rate",
    columns = contains("rejrate")
  ) |>
  tab_spanner(
    label = md("$X^2$ (df)"),
    columns = contains("X2df")
  ) |>
  cols_label(
    n = "Sample size",
    contains("PML") ~ "PML",
    contains("PMLW") ~ "PMLW"
  )

# Save tables ------------------------------------------------------------------
save(tab_bias_gt, tab_se_gt, tab_test_stats_gt, file = "sims_uneq_samp_mod5.RData")
