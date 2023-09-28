# This script to check the bias of the unequal probability sampling aka
# informative sampling.

library(tidyverse)
theme_set(theme_bw())
library(lavaan)
library(survey)
library(lavaan.bingof)
library(furrr)
library(kableExtra)
library(progressr)

# Informative sample -----------------------------------------------------------
model_no <- 1
dat <- gen_data_bin_wt(model_no = model_no, n = 1000, seed = NULL)
fit <- sem(model = txt_mod(model_no), data = dat, estimator = "PML",
           std.lv = TRUE, sampling.weights = "wt")
true_vals <- get_true_values(model_no, arrange = c("lambda", "tau", "rho"))
(bias <- lavaan:::coef(fit) - true_vals)
all_tests(fit)

make_test_stats <- function(model_no = 1, n = 1000, Sigma2_method = NULL) {
  dat <- gen_data_bin_wt(model_no, n = n, seed = NULL)
  fit <- sem(model = txt_mod(model_no), data = select(dat, -wt), estimator = "PML",
             std.lv = TRUE, sampling.weights = NULL) %>%
    suppressWarnings()

  bind_rows(
    bind_cols(all_tests(fit, Sigma2 = "theoretical"), Sigma2 = "theoretical"),
    bind_cols(all_tests(fit, Sigma2 = "weighted"), Sigma2 = "weighted"),
    bind_cols(all_tests(fit, Sigma2 = "force_unweighted"), Sigma2 = "unweighted")
  )
}

plan(multisession, workers = 30)
nsims <- 250
res <- list()
for (model_no in 1) {
  cat(paste("\nRunning model", model_no, "\n"))
  possfn <- possibly(make_test_stats, NA)

  out <-
    future_map(seq_len(nsims), \(x) {
      possfn(model_no, n = 5000)
    }, .progress = TRUE, .options = furrr_options(seed = NULL))

  out <-
    out[sapply(out, is_tibble)] %>%
    do.call(rbind, .) %>%
    mutate(model_no = model_no)

  res <- c(res, list(out))
  cat("\n")
}
res_df <- do.call(rbind, res)
res_df %>%
  group_by(model_no, name, Sigma2) %>%
  summarise(rej_rate = mean(pval < 0.05), .groups = "drop") %>%
  pivot_wider(id_cols = c(model_no, name), names_from = Sigma2,
              values_from = rej_rate) %>%
  print(n = 100)
















# Stratified sample ------------------------------------------------------------
model_no <- 1
pop <- make_population(model_no)
dat <- gen_data_bin_strat(pop, n = 1000, seed = NULL)
fit <- sem(model = txt_mod(model_no), data = dat, estimator = "PML",
           std.lv = TRUE, sampling.weights = "wt")
all_tests(fit)
all_tests(fit, Sigma2 = lavaan.bingof:::create_Sigma2_matrix_complex_x(fit))
all_tests(fit, Sigma2 = lavaan.bingof:::create_Sigma2_matrix_complex_x(fit, "strat"))

# Testing Type I error for weighted samples ------------------------------------
make_test_stats <- function(model_no = 1, n = 1000, pop) {
  if(missing(pop)) pop <- make_population(model_no, seed = NULL)
  dat <- gen_data_bin_strat(pop, n = n, seed = NULL)
  fit <- sem(model = txt_mod(model_no), data = dat, estimator = "PML",
             std.lv = TRUE, sampling.weights = "wt") %>%
    suppressWarnings()

  all_tests(fit, Sigma2 = lavaan.bingof:::create_Sigma2_matrix_complex_x(fit, "strat"))
}

plan(multisession, workers = 30)
nsims <- 100
res <- list()
for (model_no in 1:5) {
  cat(paste("\nRunning model", model_no, "\n"))
  pop <- make_population(model_no)
  possfn <- possibly(make_test_stats, NA)

  out <-
    future_map(seq_len(nsims), \(x) {
      possfn(model_no, n = 10000)
    }, .progress = TRUE, .options = furrr_options(seed = NULL))

  out <-
    out[sapply(out, is_tibble)] %>%
    do.call(rbind, .) %>%
    mutate(model_no = model_no)

  res <- c(res, list(out))
}
res_df <- do.call(rbind, res)
res_df %>%
  group_by(model_no, name) %>%
  summarise(rej_rate = mean(pval < 0.05), .groups = "drop") %>%
  pivot_wider(id_cols = name, names_from = model_no, values_from = rej_rate)




















make_bias_table <- function(n = 1000, .pop = pop, samp = "wtd") {
  samp <- match.arg(samp, c("strat", "clust", "strcl", "wtd"))

  if (samp == "wtd")   dat <- gen_data_bin_wt(model_no = model_no, n = n)
  if (samp == "strat") dat <- gen_data_bin_strat(.pop, n = n)
  if (samp == "clust") dat <- gen_data_bin_clust(.pop, n = n)
  if (samp == "strcl") dat <- gen_data_bin_strcl(.pop, n = n)

  mod <- txt_mod(model_no)

  # Ignore weights
  fit_ml <- lavaan::sem(model = mod, data = dat, std.lv = TRUE)  # DWLS
  fit_pl <- lavaan::sem(model = mod, data = dat, std.lv = TRUE,
                        estimator = "PML")
  # Fit with weights
  fit_mlw <- lavaan::sem(model = mod, data = dat, std.lv = TRUE,
                         sampling.weights = "wt")  # DWLS
  fit_plw <- lavaan::sem(model = mod, data = dat, std.lv = TRUE,
                         sampling.weights = "wt", estimator = "PML")

  true_vals <- get_true_values(model_no, arrange = c("lambda", "tau", "rho"))

  tibble(
    name = names(coef(fit_ml)),
    kind = names(true_vals) %>% gsub("\\d", "", .),
    ml  = coef(fit_ml),
    pl  = coef(fit_pl),
    mlw = coef(fit_mlw),
    plw = coef(fit_plw),
    truth = true_vals
  ) %>%
    # Get absolute bias
    mutate(across(all_of(c("ml", "pl", "mlw", "plw")), \(x) sqrt((x - truth) ^ 2)))
}


new_levels <- c(paste0("eta2=~y", 6:10), paste0("eta3=~y", 11:15))
names(new_levels) <- c(paste0("eta1=~y", 6:10), paste0("eta1=~y", 11:15))

# Summarise in a table
options(knitr.kable.NA = "")
res_df %>%
  group_by(model_no, name, kind) %>%
  summarise(across(c(pl, plw), mean), .groups = "drop") %>%
  # summarise(across(c(pl, plw), \(x) paste0(iprior::dec_plac(mean(x), 2),
  #                                          " (",
  #                                          iprior::dec_plac(sd(x), 2),
  #                                          ")"))) %>%
  mutate(name = factor(name, levels = c(paste0("eta1=~y", 1:15),
                                        paste0("eta2=~y", 6:10),
                                        paste0("eta3=~y", 11:15),
                                        paste0("y", 1:15, "|t1"),
                                        "eta1~~eta2", "eta1~~eta3", "eta2~~eta3"
  ))) %>%
  mutate(name = fct_recode(name, !!!new_levels)) %>%
  pivot_wider(id_cols = c(name, kind), names_from = model_no,
              values_from = c(pl, plw)) %>%
  # filter(kind != "lambda") %>%
  mutate(kind = case_when(kind == "lambda" ~ "Loadings",
                          kind == "tau" ~ "Thresholds",
                          kind == "rho" ~ "Fac. corr.")) %>%
  select(kind, name,
         pl_1, plw_1,
         pl_2, plw_2,
         pl_3, plw_3,
         pl_4, plw_4,
         pl_5, plw_5) %>%
  arrange(name) %>%
  kbl(format = "latex", digits = 2, booktabs = TRUE)  %>%
  add_header_above(c(" " = 2, "Model 1" = 2, "Model 2" = 2, "Model 3" = 2,
                     "Model 4" = 2, "Model 5" = 2)) %>%
  collapse_rows(1:2, latex_hline = "none", row_group_label_position = 'stack')

res_df %>%
  filter(model_no == 3) %>%
  pivot_longer(cols = c("pl", "plw"),
               names_to = "type", values_to = "bias") %>%
  mutate(type = factor(type, levels = c("pl", "plw")),
         kind = factor(kind, levels = c("lambda", "tau", "rho"))) %>%
  ggplot(aes(name, bias, fill = type)) +
  # geom_point(position = position_dodge(width = 0.5)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid( ~ kind, scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Checking bias of estimated probabilities -------------------------------------
compare_model_probs <- function(.model_no = 1, n = 1000, samp = "wtd",
                                .seed = NULL) {

  if (samp == "wtd") {
    dat <- gen_data_bin_wt(model_no = .model_no, n = n, seed = .seed)
    fit <- sem(model = txt_mod(.model_no), data = dat, estimator = "PML",
               std.lv = TRUE, sampling.weights = "wt") %>%
      suppressWarnings()
    svy <- svydesign(ids = ~ 1, weights = ~ wt, data = dat)
  }
  if (samp == "strat") {
    dat <- gen_data_bin_strat(make_population(.model_no, seed = NULL), n = n, seed = .seed)
    fit <- sem(model = txt_mod(.model_no), data = dat, estimator = "PML",
               std.lv = TRUE, sampling.weights = "wt") %>%
      suppressWarnings()
    svy <- svydesign(ids = ~ 0, strata = ~ type, weights = ~ wt, data = dat)
  }
  if (samp == "clust") {
    dat <- gen_data_bin_clust(make_population(.model_no, seed = NULL), n = n, seed = .seed)
    fit <- sem(model = txt_mod(.model_no), data = dat, estimator = "PML",
               std.lv = TRUE, sampling.weights = "wt") %>%
      suppressWarnings()
    svy <- svydesign(ids = ~ school + class, weights = ~ wt, data = dat)
  }
  if (samp == "strcl") {
    dat <- gen_data_bin_strcl(make_population(.model_no, seed = NULL), n = n, seed = .seed)
    fit <- sem(model = txt_mod(.model_no), data = dat, estimator = "PML",
               std.lv = TRUE, sampling.weights = "wt") %>%
      suppressWarnings()
    svy <- svydesign(ids = ~ school + class, strata = ~ type,
                     weights = ~ wt, data = dat, nest = TRUE)
  }
  if (samp == "srs") {
    dat <- gen_data_bin(model_no = .model_no, n = n, seed = .seed)
    fit <- sem(model = txt_mod(.model_no), data = dat, estimator = "PML",
               std.lv = TRUE) %>%
      suppressWarnings()
    svy <- NULL
  }

  Omega2 <- lavaan.bingof:::calc_test_stuff(fit, svy)$Omega2
  tmp <- eigen(Omega2)
  Omegahalfinv <- with(tmp, vectors %*% diag(1 / sqrt(values)) %*% t(vectors))

  truth <- with(lavaan.bingof:::get_theoretical_uni_bi_moments(.model_no),
                c(pidot1, pidot2))
  res <- list(truth = truth) %>%
    c(., with(lavaan.bingof:::get_uni_bi_moments(fit),
              list(wtdprop = c(pdot1, pdot2),
                   pihat = c(pidot1, pidot2))),
      with(lavaan.bingof:::get_uni_bi_moments(fit, wtd = FALSE),
           list(prop = c(pdot1, pdot2))))

  p <- nrow(lavaan.bingof:::get_Lambda(.model_no))
  ynames <- c(
    paste0("y", 1:p),
    as_tibble(t(combn(p, 2))) %>%
      rename(i = 1, j = 2) %>%
      unite("ij", sep = ",") %>%
      pull() %>%
      paste0("y", .)
  )

  bind_cols(tibble(idx = ynames), (res)) %>%
    mutate(idx = factor(idx, levels = idx)) %>%
    mutate(g2t = c(Omegahalfinv %*% (wtdprop - truth)) / sqrt(n)) %>%
    mutate(g2 = c(Omegahalfinv %*% (wtdprop - pihat)) / sqrt(n)) %>%
    mutate(seed = .seed)
}; compare_model_probs(n = 1000, samp = "strat")

res <- list()
for (samp_method in c( "strat", "clust", "strcl")) {
  for (model_no in 1) {
    for (samp_size in c(1000L, 5000L, 100000L)) {
      plan(multisession, workers = 25)
      cat(paste0("\nMethod: ", samp_method, " / Model no: ", model_no, " / Sample size: ", samp_size, "\n"))

      posscmp <- possibly(compare_model_probs, tibble(idx = NA, truth = NA,
                                                      wtdprop = NA, pihat = NA,
                                                      prop = NA, seed = NA))

      out <-
        future_map(
          1:50, ~ posscmp(model_no, samp_size, samp_method, .x * 2),
          .progress = TRUE, .options = furrr_options(seed = NULL)
        ) %>%
        do.call(rbind, .) %>%
        mutate(model_no = model_no, n = paste0("n = ", as.character(samp_size)),
               samp = samp_method)

      res <- c(res, list(out))
    }
  }
}

res_df <- do.call(rbind, res) %>%
  # mutate(n = fct_recode(n, "n = 100000" = "n = 1e+05")) %>%
  mutate(n = fct_relevel(n, "n = 100000", after = Inf))

res_df %>%
  pivot_longer(c(pihat, wtdprop, prop), names_to = "est") %>%
  unite(linegrp, seed, est, remove = FALSE) %>%
  mutate(samp = fct_relevel(samp, c("srs", "wtd", "clust"))) %>%
  mutate(samp = fct_recode(samp, "Independent sample (SRS)" = "srs",
                           "Weighted sample" = "wtd",
                           "Cluster sample" = "clust")) %>%
  mutate(est = fct_relevel(est, c("prop", "wtdprop", "pihat"))) %>%
  mutate(est = fct_recode(est, `1` = "prop", `2` = "wtdprop", `3` = "pihat")) %>%
  ggplot() +
  geom_line(aes(as.numeric(idx), value, group = linegrp, col = est),
            alpha = 0.5, linewidth = 0.2) +
  geom_line(data = res_df %>% group_by(idx) %>% summarise(value = mean(truth)),
            aes(as.numeric(idx), value, col = "4")) +
  # geom_point(position = position_dodge(width = 0.2)) +
  facet_grid(n ~ samp) +
  guides(col = guide_legend(override.aes = list(linewidth = 1)))  +
  scale_colour_manual(lab = c("Unwtd. prop.", "Wtd. prop", "Model prob.", "Truth"),
                      values = c(iprior::gg_col_hue(3), "black")) +
  scale_x_continuous(breaks = 1:15, labels = unique(res_df$idx)) +
  labs(x = NULL, y = "Probability", col = NULL) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

res_df %>%
  pivot_longer(c(pihat, wtdprop, prop), names_to = "est") %>%
  unite(linegrp, seed, est, remove = FALSE) %>%
  mutate(samp = fct_relevel(samp, c("srs", "wtd", "clust"))) %>%
  mutate(samp = fct_recode(samp, "Independent sample (SRS)" = "srs",
                           "Weighted sample" = "wtd",
                           "Cluster sample" = "clust")) %>%
  mutate(est = fct_relevel(est, c("prop", "wtdprop", "pihat"))) %>%
  mutate(est = fct_recode(est, `1` = "prop", `2` = "wtdprop", `3` = "pihat")) %>%
  mutate(bias = value - truth) %>%
  ggplot() +
  geom_line(aes(as.numeric(idx), bias, group = linegrp, col = est),
            alpha = 0.5, linewidth = 0.2) +
  facet_grid(n ~ samp) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  guides(col = guide_legend(override.aes = list(linewidth = 1))) +
  scale_colour_manual(lab = c("Unwtd. prop.", "Wtd. prop", "Model prob.", "Truth"),
                      values = c(iprior::gg_col_hue(3))) +
  scale_x_continuous(breaks = 1:15, labels = unique(res_df$idx)) +
  labs(x = NULL, y = "Bias", col = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

res_df %>%
  pivot_longer(c(g2), names_to = "est") %>%
  unite(linegrp, seed, est, remove = FALSE) %>%
  mutate(samp = fct_relevel(samp, c("srs", "wtd", "clust"))) %>%
  mutate(samp = fct_recode(samp, "Independent sample (SRS)" = "srs",
                           "Weighted sample" = "wtd",
                           "Cluster sample" = "clust")) %>%
  # mutate(est = fct_relevel(est, c("prop", "wtdprop", "pihat"))) %>%
  # mutate(est = fct_recode(est, `1` = "prop", `2` = "wtdprop", `3` = "pihat")) %>%
  ggplot() +
  geom_line(aes(as.numeric(idx), value, group = linegrp, col = est),
            alpha = 0.5, linewidth = 0.2) +
  facet_grid(n ~ samp) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  guides(col = guide_legend(override.aes = list(linewidth = 1))) +
  # scale_colour_manual(lab = c("Unwtd. prop.", "Wtd. prop", "Model prob.", "Truth"),
  #                     values = c(iprior::gg_col_hue(3))) +
  scale_x_continuous(breaks = 1:15, labels = unique(res_df$idx)) +
  labs(x = NULL, y = "Value", col = NULL) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))























