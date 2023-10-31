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

plan(multisession, workers = 8)
nsims <- 100

calc_X2 <- function(fit) {
  tab <- lavaan::lavTables(fit, dimension = 0L) |>
    suppressWarnings()
  N <- tab$nobs[1]
  X <- as.data.frame(fit@Data@X[[1]])
  wt <- fit@Data@weights[[1]]

  # Calc weighted props
  if (!is.null(wt)) {
    # N <- lavaan.bingof:::neff(wt)
    wtprops <-
      as_tibble(X) |>
      unite("pattern", everything(), remove = TRUE, sep = "") |>
      mutate(wt = wt) |>
      group_by(pattern) |>
      summarise(wtfreq = sum(wt)) |>
      mutate(wtprop = wtfreq / sum(wtfreq))

    X2 <- with(left_join(tab, wtprops, by = "pattern"), {
      N * sum((est.prop - wtprop) ^ 2 / est.prop)
    })
  } else {
    X2 <- with(tab, {
      N * sum((est.prop - obs.prop) ^ 2 / est.prop)
    })
  }

  # p <- 5
  # q <- 1
  # # df <- 2 ^ p - 1 - p + p * q - q * (q - 1) / 2
  df <- 2^5 - 1 - 10

  list(X2 = X2, df = df)
}

make_wt_res <- function(samp_size = 5000, type = "bias_se", .popfac = 1 / 0.01) {
  model_no <- 1
  dat <- gen_data_bin_wt(model_no, n = samp_size, popfac = .popfac)
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
                            est    = se0$est,
                            se     = se0$se,
                            se2    = se)),
      bind_cols(tmp, tibble(method = "PMLW",
                            est    = se1$est,
                            se     = se1$se,
                            se2    = se * sqrt(samp_size / ntilde)
                            ))
    ) %>%
      mutate(
        unc = qnorm(0.975) * se,
        covered = (truth >= est - unc) & (truth <= est + unc),
        covered2 = (truth >= est - qnorm(0.975) * se2) &
          (truth <= est + qnorm(0.975) * se2),
      )
  }
  if (type == "test_stats") {
    # Test statistics ----------------------------------------------------------
    # PML
    tmp <- lavaan::fitMeasures(fit0)
    tab0 <-
      Pearson_test(fit0, Sigma2 = "theoretical") %>%
      mutate(chisq = calc_X2(fit0)$X2, #tmp["chisq"],
             df2 =  calc_X2(fit0)$df, #tmp["df"],
             pval2 = pchisq(chisq, df, lower.tail = FALSE))

    tmp <- lavaan::fitMeasures(fit1)
    tab1 <-
      Pearson_test(fit1, Sigma2 = "theoretical") %>%
      mutate(chisq = calc_X2(fit1)$X2, #tmp["chisq"],
             df2 =  calc_X2(fit1)$df, #tmp["df"],
             pval2 = pchisq(chisq, df, lower.tail = FALSE))

    res <- bind_rows(
      bind_cols(method = "PML", tab0),
      bind_cols(method = "PMLW", tab1),
    )
  }

  as_tibble(res)

}
possfn <- possibly(make_wt_res, NA)

# Bias and se sims -------------------------------------------------------------

res <-
  future_map(seq_len(nsims), \(x) {
    possfn(type = "bias_se")
  }, .progress = TRUE, .options = furrr_options(seed = NULL))

res_df_biasse <-
  res[sapply(res, is_tibble)] %>%
  do.call(rbind, .)

# Bias table
# res_df %>%
#   mutate(bias = est - truth) %>%
#   group_by(method, name, kind, truth) %>%
#   summarise(bias = mean((bias)), .groups = "drop") %>%
#   pivot_wider(id_cols = c(name, kind, truth), names_from = method,
#               values_from = bias) %>%
#   kbl(format = "latex", digits = 2, booktabs = TRUE)

res_df_biasse %>%
  ggplot(aes(est, col = method, group = method)) +
  geom_density() +
  facet_wrap(name ~ ., scales = "free")

# SE table
res_df_biasse %>%
  mutate(bias = est - truth) %>%
  group_by(kind, name, method, truth) %>%
  summarise(bias = mean(bias),
            cov = mean(covered),
            cov2 = mean(covered2),
            # sdse = sd(est) / mean(se),
            sd = sqrt(mean((est - truth) ^ 2)),
            # sd = sd(est),
            se = mean(se),
            se2 = mean(se2),
            .groups = "drop") %>%
  mutate(sdse = sd / se,
         sdse2 = sd / se2) %>%
  pivot_wider(id_cols = c(kind, name, truth), names_from = method,
              values_from = c(bias, cov, sdse)) %>%
  mutate(kind = rep(c("Loadings", "Thresholds"), each = 5),
         name = c(paste0("$\\lambda_", 1:5, "$"),
                  paste0("$\\tau_", 1:5, "$"))) %>%
  kbl(format = "latex", digits = 2, booktabs = TRUE, escape = FALSE) %>%
  add_header_above(c(" " = 3, "Bias" = 2, "Coverage" = 2, "SD/SE" = 2)) %>%
  collapse_rows(1:2, row_group_label_position = "stack", latex_hline = "none")
















# Test statistics sims ---------------------------------------------------------

res <- list()
for (n in c(500, 1000, 2500, 5000, 10000)) {  #
  cat(paste("\nRunning with sample size n =", n, "\n"))

  out <-
    future_map(seq_len(nsims), \(x) {
      possfn(samp_size = n, type = "test_stats", .popfac = 500000 / n)
    }, .progress = TRUE, .options = furrr_options(seed = NULL))

  out <-
    out[sapply(out, is_tibble)] %>%
    do.call(rbind, .) %>%
    mutate(n = n)

  res <- c(res, list(out))
  cat("\n")
}
res_df <- do.call(rbind, res)

res_df %>%
  group_by(n, method) %>%
  summarise(rej_rate = mean(pval < 0.05),
            X2 = iprior::dec_plac(mean(X2), 2),
            df = iprior::dec_plac(mean(df), 2),
            rej_rate2 = mean(pval2 < 0.05),
            chisq = iprior::dec_plac(mean(chisq), 2),
            df2 = iprior::dec_plac(mean(df2), 2),
            .groups = "drop") %>%
  mutate(X2df = paste0(X2, " (", df, ")"),
         chi2df = paste0(chisq, " (", df2, ")")) %>%
  pivot_wider(id_cols = n, names_from = method,
              # values_from = c(rej_rate2, chi2df)) %>%
              values_from = c(rej_rate, X2df))

  kbl(format = "latex", digits = 2, booktabs = TRUE) %>%
  add_header_above(c(" " = 1, "Rejection rate" = 2, "X2 (df)" = 2))
