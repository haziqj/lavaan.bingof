library(tidyverse)
library(lavaan.bingof)
library(furrr)
library(cli)
library(gt)
plan(multisession, workers = parallel::detectCores() - 2)
nsims <- 1000

get_clust_se <- function(object, dat) {
  
  # Obtain scores
  Hinv <- lavaan.bingof:::get_sensitivity_inv_mat(object)  # H^{-1}
  Delta <- lavaan:::computeDelta(lavmodel = object@Model)[[1]]
  SC <- lavaan:::pml_deriv1(
    Sigma.hat = object@implied$cov[[1]],
    Mu.hat = object@implied$mean[[1]],
    TH = object@implied$th[[1]],
    th.idx = object@Model@th.idx[[1]],
    num.idx = object@Model@num.idx[[1]],
    X = object@Data@X[[1]],
    eXo = NULL,
    wt = NULL,
    PI = NULL,
    lavcache = object@Cache[[1]],
    missing = object@Data@missing,
    scores = TRUE,
    negative = FALSE
  )
  
  # Combine with weights information from dat
  SCwt <-
    dat |>
    select(school, wt) |>
    bind_cols(as_tibble(SC))
  
  # Get the clustered J matrix and Godambe matrix
  clusters <- unique(SCwt$school)
  nclust <- length(clusters)
  zb <- list()
  for (b in seq_along(clusters)) {
    dat_b <- 
      SCwt |>
      filter(school == clusters[b])
    
    tmp <- select(dat_b, -school, -wt) * dat_b$wt
    tmp <- apply(tmp, 2, sum)
    zb[[b]] <- tmp
  }
  zbar <- apply(do.call(cbind, zb), 1, mean)
  
  B1c <-
    lapply(zb, \(z) tcrossprod(z - zbar)) |>
    Reduce(f = `+`)
  Jc <- nclust / (nclust - 1) * (t(Delta) %*% B1c %*% Delta) / nrow(dat)
  Gc <- Hinv %*% Jc %*% Hinv  # inverse Godambe matrix
  
  # Return standard errors
  out <- sqrt(diag(Gc) / nrow(dat)) 
  names(out) <- names(coef(object))
  out
}

make_clust_res <- function(samp_size = 1000, type = "bias_se") {
  model_no <- 1
  pop <- make_population(model_no, seed = 123)
  dat <- gen_data_bin_clust(pop, n = samp_size)
  fit <- lavaan::sem(model = txt_mod(model_no), data = dat, estimator = "PML",
                     std.lv = TRUE, sampling.weights = "wt")
  
  if (type == "bias_se") {
    # Bias and se results ------------------------------------------------------
    true_vals <- get_true_values(model_no, arrange = c("lambda", "tau", "rho"))
    se0 <- se1 <- lavaan::parTable(fit) |> filter(free > 0)
    se1$se <- get_clust_se(fit, dat)
    
    tab <- tibble(
      name  = names(lavaan::coef(fit)),
      kind  = gsub("\\d", "", names(true_vals)),
      truth = true_vals,
      n     = samp_size
    )
    res <- 
      bind_rows(
        bind_cols(tab, tibble(method = "PMLW",
                              est    = se0$est,
                              se     = se0$se)),
        bind_cols(tab, tibble(method = "PMLW (Clust)",
                              est    = se1$est,
                              se     = se1$se))
      ) |>
      mutate(
        unc = qnorm(0.975) * se,
        covered = (truth >= est - unc) & (truth <= est + unc)
      )
  }

  # Check for Heywood cases
  res$var_ok <- lavaan:::lav_object_post_check(fit)
  res$converged <- fit@Fit@converged
  as_tibble(res)
}
make_clust_res <- possibly(make_clust_res, NA)

# Run simulations --------------------------------------------------------------

# Bias, coverage and SE
res <- list()
for (n in c(500, 1000, 5000)) {
  cli_alert_info("Running with sample size n = {n}")
  
  out <-
    future_map(seq_len(nsims), \(x) {
      make_clust_res(samp_size = n, type = "bias_se")
    }, .progress = TRUE, .options = furrr_options(seed = NULL))
  
  out <-
    out[sapply(out, is_tibble)] |>
    do.call(what = rbind)
  
  res <- c(res, list(out))
  cat("\n")
}
res_bias_se_df <- do.call(rbind, res)


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
  mutate(
    name = gsub("eta1=~y", "$\\\\lambda_", name),
    name = gsub("y", "$\\\\tau_", name),
    name = gsub("\\|t1", "", name),
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
    label = md("**Thresholds**"),
    rows = contains("tau")
  ) |>  
  tab_row_group(
    label = md("**Loadings**"),
    rows = contains("lambda")
  ) |>
  cols_width(truth ~ px(70)) |>
  cols_label(
    contains("PML_") ~ "PML",
    contains("PMLW_") ~ "PMLW",
    truth = "True values"
  )

# Create se table --------------------------------------------------------------
tab_clust_se <-
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

tab_clust_se_gt <-
  tab_clust_se |>
  mutate(
    name = gsub("eta1=~y", "$\\\\lambda_", name),
    name = gsub("y", "$\\\\tau_", name),
    name = gsub("\\|t1", "", name),
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
    label = md("**Thresholds**"),
    rows = contains("tau")
  ) |>  
  tab_row_group(
    label = md("**Loadings**"),
    rows = contains("lambda")
  ) |>
  cols_label(
    TRUE ~ "PMLW (Clust)",
    contains("PMLW_") ~ "PMLW"
  )

# Save tables ------------------------------------------------------------------
save(tab_clust_se_gt, file = "simulations/sims_clust.RData")
