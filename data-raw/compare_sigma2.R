library(tidyverse)
library(lavaan.bingof)
library(furrr)
plan(multisession, workers = parallel::detectCores() - 2)

# For each of the samples obtained using 1) wtd, 2) srs, 3) strat, 4) clust, and
# 5) strcl, we will compare the performance of the LIGOF tests for Model 5
# specifically. Note that each of these sampling we can either use weights or
# ignore the weights (except for srs). Then, we compare the use of these
# estimators of Sigma2:
#
# a. theoretical
# b. force_unweighted
# c. weighted
# d. strat
# e. population

model_no <- 3
samp_size <- 1000
no_sims <- 1000
alt_hyp <- FALSE
pop <- make_population(model_no, seed = 21324, Sigma2_attr = TRUE, H1 = alt_hyp)

run_compare_sigma2_sims <- function(
    samp = c("wtd", "srs", "strat", "clust", "strcl"),
    B = 10) {
  samp <- match.arg(samp, c("wtd", "srs", "strat", "clust", "strcl"))
  set.seed(NULL)

  if (samp == "wtd") {
    the_fun <- possibly(sim_wtd, NA)
  }
  if (samp == "srs") {
    the_fun <- possibly(sim_srs, NA)
  }
  if (samp == "strat") {
    the_fun <- possibly(sim_strat, NA)
  }
  if (samp == "clust") {
    the_fun <- possibly(sim_clust, NA)
  }
  if (samp == "strcl") {
    the_fun <- possibly(sim_strcl, NA)
  }

  cli::cli_alert_info("Now running: {samp}")

  out <- furrr::future_map(
    .x = 1:B,
    .f = the_fun,
    .options = furrr::furrr_options(seed = TRUE),
    .progress = TRUE
  )
  cat("\n")
  out
}

sim_wtd <- function(i) {
  dat <- gen_data_bin_wt(model_no, n = samp_size, H1 = alt_hyp)

  # Unweighted results ---------------------------------------------------------
  fit0 <- lavaan::sem(
    model = txt_mod(model_no), data = dat, estimator = "PML",
    std.lv = TRUE
  )

  res_nowt <-
    tribble(
      ~sigma2, ~res,
      "theoretical", all_tests(fit0, sim = i, Sigma2 = "theoretical"),
      "multinomial", all_tests(fit0, sim = i, Sigma2 = "multinomial"),
      # "weighted"        , NA,
      # "strat"           , NA,
      "force_unweighted", all_tests(fit0, sim = i, Sigma2 = "force_unweighted")
    )

  # Weighted results -----------------------------------------------------------
  fit1 <- lavaan::sem(
    model = txt_mod(model_no), data = dat, estimator = "PML",
    std.lv = TRUE, sampling.weights = "wt"
  )

  res_wt <-
    tribble(
      ~sigma2, ~res,
      "theoretical", all_tests(fit1, sim = i, Sigma2 = "theoretical"),
      "multinomial", all_tests(fit1, sim = i, Sigma2 = "multinomial"),
      "weighted", all_tests(fit1, sim = i, Sigma2 = "weighted"),
      # "strat"           , NA,
      "force_unweighted", all_tests(fit1, sim = i, Sigma2 = "force_unweighted")
    )

  # Return results -------------------------------------------------------------
  bind_rows(
    bind_cols(wt = "wt", res_wt),
    bind_cols(wt = "nowt", res_nowt)
  )
}

sim_srs <- function(i) {
  dat <- gen_data_bin_srs(pop, n = samp_size)

  # Unweighted results ---------------------------------------------------------
  fit0 <- lavaan::sem(
    model = txt_mod(model_no), data = dat, estimator = "PML",
    std.lv = TRUE
  )
  res_nowt <-
    tribble(
      ~sigma2, ~res,
      "theoretical", all_tests(fit0, sim = i, Sigma2 = "theoretical"),
      "multinomial", all_tests(fit0, sim = i, Sigma2 = "multinomial"),
      # "weighted"        , all_tests(fit0, sim = i, Sigma2 = "weighted"),
      # "strat"           , all_tests(fit0, sim = i, Sigma2 = "strat"),
      "force_unweighted", all_tests(fit0, sim = i, Sigma2 = "force_unweighted"),
      "population", all_tests(fit0, sim = i, Sigma2 = attr(pop, "Sigma2"))
    )

  # Weighted results -----------------------------------------------------------
  # fit1 <- lavaan::sem(model = txt_mod(model_no), data = dat, estimator = "PML",
  #                     std.lv = TRUE, sampling.weights = "wt")

  res_wt <-
    tribble(
      ~sigma2, ~res,
      # "theoretical"     , all_tests(fit1, sim = i, Sigma2 = "theoretical"),
      # "weighted"        , all_tests(fit1, sim = i, Sigma2 = "weighted"),
      # "strat"           , all_tests(fit1, sim = i, Sigma2 = "strat"),
      # "force_unweighted", all_tests(fit1, sim = i, Sigma2 = "force_unweighted"),
      # "population"      , all_tests(fit1, sim = i, Sigma2 = attr(pop, "Sigma2"))
    )

  # Return results -------------------------------------------------------------
  bind_rows(
    bind_cols(wt = "wt", res_wt),
    bind_cols(wt = "nowt", res_nowt)
  )
}

sim_strat <- function(i) {
  dat <- gen_data_bin_strat(pop, n = samp_size)

  # Unweighted results ---------------------------------------------------------
  fit0 <- lavaan::sem(
    model = txt_mod(model_no), data = dat, estimator = "PML",
    std.lv = TRUE
  )
  res_nowt <-
    tribble(
      ~sigma2, ~res,
      "theoretical", all_tests(fit0, sim = i, Sigma2 = "theoretical"),
      "multinomial", all_tests(fit0, sim = i, Sigma2 = "multinomial"),
      # "weighted"        , all_tests(fit0, sim = i, Sigma2 = "weighted"),
      # "strat"           , all_tests(fit0, sim = i, Sigma2 = "strat"),
      "force_unweighted", all_tests(fit0, sim = i, Sigma2 = "force_unweighted"),
      "population", all_tests(fit0, sim = i, Sigma2 = attr(pop, "Sigma2"))
    )

  # Weighted results -----------------------------------------------------------
  fit1 <- lavaan::sem(
    model = txt_mod(model_no), data = dat, estimator = "PML",
    std.lv = TRUE, sampling.weights = "wt"
  )

  res_wt <-
    tribble(
      ~sigma2, ~res,
      "theoretical", all_tests(fit1, sim = i, Sigma2 = "theoretical"),
      "multinomial", all_tests(fit1, sim = i, Sigma2 = "multinomial"),
      "weighted", all_tests(fit1, sim = i, Sigma2 = "weighted"),
      "strat", all_tests(fit1, sim = i, Sigma2 = "strat"),
      "force_unweighted", all_tests(fit1, sim = i, Sigma2 = "force_unweighted"),
      "population", all_tests(fit1, sim = i, Sigma2 = attr(pop, "Sigma2"))
    )

  # Return results -------------------------------------------------------------
  bind_rows(
    bind_cols(wt = "wt", res_wt),
    bind_cols(wt = "nowt", res_nowt)
  )
}

sim_clust <- function(i) {
  dat <- gen_data_bin_clust(pop, n = samp_size)

  # Unweighted results ---------------------------------------------------------
  fit0 <- lavaan::sem(
    model = txt_mod(model_no), data = dat, estimator = "PML",
    std.lv = TRUE
  )
  res_nowt <-
    tribble(
      ~sigma2, ~res,
      "theoretical", all_tests(fit0, sim = i, Sigma2 = "theoretical"),
      "multinomial", all_tests(fit0, sim = i, Sigma2 = "multinomial"),
      # "weighted"        , all_tests(fit0, sim = i, Sigma2 = "weighted"),
      # "strat"           , all_tests(fit0, sim = i, Sigma2 = "strat"),
      "force_unweighted", all_tests(fit0, sim = i, Sigma2 = "force_unweighted"),
      "population", all_tests(fit0, sim = i, Sigma2 = attr(pop, "Sigma2"))
    )

  # Weighted results -----------------------------------------------------------
  fit1 <- lavaan::sem(
    model = txt_mod(model_no), data = dat, estimator = "PML",
    std.lv = TRUE, sampling.weights = "wt"
  )

  res_wt <-
    tribble(
      ~sigma2, ~res,
      "theoretical", all_tests(fit1, sim = i, Sigma2 = "theoretical"),
      "multinomial", all_tests(fit1, sim = i, Sigma2 = "multinomial"),
      "weighted", all_tests(fit1, sim = i, Sigma2 = "weighted"),
      # "strat"           , all_tests(fit1, sim = i, Sigma2 = "strat"),
      "force_unweighted", all_tests(fit1, sim = i, Sigma2 = "force_unweighted"),
      "population", all_tests(fit1, sim = i, Sigma2 = attr(pop, "Sigma2"))
    )

  # Return results -------------------------------------------------------------
  bind_rows(
    bind_cols(wt = "wt", res_wt),
    bind_cols(wt = "nowt", res_nowt)
  )
}

sim_strcl <- function(i) {
  dat <- gen_data_bin_strcl(pop, n = samp_size)

  # Unweighted results ---------------------------------------------------------
  fit0 <- lavaan::sem(
    model = txt_mod(model_no), data = dat, estimator = "PML",
    std.lv = TRUE
  )
  res_nowt <-
    tribble(
      ~sigma2, ~res,
      "theoretical", all_tests(fit0, sim = i, Sigma2 = "theoretical"),
      "multinomial", all_tests(fit0, sim = i, Sigma2 = "multinomial"),
      # "weighted"        , all_tests(fit0, sim = i, Sigma2 = "weighted"),
      # "strat"           , all_tests(fit0, sim = i, Sigma2 = "strat"),
      "force_unweighted", all_tests(fit0, sim = i, Sigma2 = "force_unweighted"),
      "population", all_tests(fit0, sim = i, Sigma2 = attr(pop, "Sigma2"))
    )

  # Weighted results -----------------------------------------------------------
  fit1 <- lavaan::sem(
    model = txt_mod(model_no), data = dat, estimator = "PML",
    std.lv = TRUE, sampling.weights = "wt"
  )

  res_wt <-
    tribble(
      ~sigma2, ~res,
      "theoretical", all_tests(fit1, sim = i, Sigma2 = "theoretical"),
      "multinomial", all_tests(fit1, sim = i, Sigma2 = "multinomial"),
      "weighted", all_tests(fit1, sim = i, Sigma2 = "weighted"),
      # "strat"           , all_tests(fit1, sim = i, Sigma2 = "strat"),
      "force_unweighted", all_tests(fit1, sim = i, Sigma2 = "force_unweighted"),
      "population", all_tests(fit1, sim = i, Sigma2 = attr(pop, "Sigma2"))
    )

  # Return results -------------------------------------------------------------
  bind_rows(
    bind_cols(wt = "wt", res_wt),
    bind_cols(wt = "nowt", res_nowt)
  )
}

res <-
  # tibble(samp = c("wtd", "srs", "strat", "clust", "strcl")) |>
  tibble(samp = c("strcl")) |>
  mutate(res = purrr::map(
    .x = samp,
    .f = \(x) run_compare_sigma2_sims(x, B = no_sims)
  ))
# save(res, file = "compare_sigma2.RData")

# res |>
#   unnest(res) |>
#   unnest(res) |>
#   unnest(res) |>
#   filter(
#     samp == "strcl",
#     # name %in% c("Wald", "Pearson,MM3")
#   ) |>
#   ggplot(aes(X2, y = name, fill = wt)) +
#   geom_violin(draw_quantiles = c(0.5)) +
#   facet_grid(. ~ sigma2)

res |>
  unnest(res) |>
  unnest(res) |>
  unnest(res) |>
  mutate(
    rej = pval < 0.05,
    ok = converged & Omega2_rank >= S
  ) |>
  summarise(
    rej_rate = mean(rej[ok]),
    .by = c(samp, wt, sigma2, name)
  ) |>
  mutate(
    sigma2 = case_when(
      sigma2 == "population" ~ "Sigma2 = Pop.",
      sigma2 == "theoretical" ~ "Sigma2 = Theor.",
      sigma2 == "multinomial" ~ "Sigma2 = Multn.",
      sigma2 == "weighted" ~ "Sigma2 = Wtd.",
      sigma2 == "strat" ~ "Sigma2 = Strat.",
      sigma2 == "force_unweighted" ~ "Sigma2 = Unwtd."
    ),
    sigma2 = factor(sigma2, levels = c("Sigma2 = Theor.", "Sigma2 = Multn.", "Sigma2 = Unwtd.", "Sigma2 = Wtd.", "Sigma2 = Strat.", "Sigma2 = Pop.")),
    wt = case_when(
      wt == "wt" ~ "Weighted",
      wt == "nowt" ~ "Ignore wt."
    ),
    samp = case_when(
      samp == "wtd" ~ "Informative",
      samp == "srs" ~ "SRS",
      samp == "strat" ~ "Stratified",
      samp == "clust" ~ "Cluster",
      samp == "strcl" ~ "Strat. + Clust."
    ),
    samp = factor(samp, levels = c("Informative", "SRS", "Stratified", "Cluster", "Strat. + Clust.")),
    name = factor(name, levels = c("Wald", "WaldVCF", "WaldDiag,MM3", "Pearson,MM3", "RSS,MM3", "Multn,MM3"))
  ) |>
  filter(sigma2 != "Sigma2 = Multn.") |>
  ggplot(aes(name, rej_rate, fill = name)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_hline(yintercept = 0.05, linetype = "dashed", col = "black", linewidth = 0.8) +
  facet_grid(samp * wt ~ sigma2) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, 0.2)) +
  labs(
    x = NULL,
    y = "Rejection rate",
    title = glue::glue("Scenario 3F15V (n={samp_size}): Rejection rates by simulation conditions"),
    caption = glue::glue("Total replications: {no_sims}")
  ) +
  guides(fill = "none") +
  ggsci::scale_fill_jama()
