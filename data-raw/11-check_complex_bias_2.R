library(tidyverse)
theme_set(theme_bw())
library(lavaan)
library(survey)
library(lavaan.bingof)
library(furrr)

model_no <- 1
pop <- make_population(model_no)

make_bias_table <- function(n = 3000, .pop = pop, samp) {
  samp <- match.arg(samp, c("strat", "clust", "strcl"))
  if (samp == "strat") dat <- gen_data_bin_strat(.pop, n = n)
  if (samp == "clust") dat <- gen_data_bin_clust(.pop, n = n)
  if (samp == "strcl") dat <- gen_data_bin_strcl(.pop, n = n)

  mod <- txt_mod(model_no)

  fit_ml <- lavaan::sem(model = mod, data = dat, std.lv = TRUE)  # DWLS
  fit_pl <- lavaan::sem(model = mod, data = dat, std.lv = TRUE,
                        estimator = "PML")
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
  )
}

res <- list()
for (model_no in 1:5) {
  for (samp in c("strat", "clust", "strcl")) {
    pop <- make_population(model_no)
    plan(multisession, workers = 30)

    out <-
      future_map(1:1000, ~make_bias_table(samp = samp), .progress = TRUE,
                 .options = furrr_options(seed = NULL)) %>%
      do.call(rbind, .) %>%
      mutate(model_no = model_no, samp = samp)

    res <- c(res, list(out))
  }
}
res <- do.call(rbind, res)
save(res, file = "vignettes/articles/check_complex_bias2.RData")


stuff <- c("ml", "pl", "mlw", "plw")
res %>%
  filter(model_no == 4) %>%
  mutate(across(all_of(stuff), \(x) abs(x - truth) ^ 1)) %>%
  pivot_longer(cols = stuff, names_to = "type",
               values_to = "bias") %>%
  mutate(type = factor(type, levels = stuff),
         kind = factor(kind, levels = c("lambda", "tau", "rho"))) %>%
  ggplot(aes(name, bias, fill = type)) +
  # geom_point(position = position_dodge(width = 0.5)) +
  geom_boxplot(outlier.shape = NA) +
  facet_grid(samp ~ ., scales = "free") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))








