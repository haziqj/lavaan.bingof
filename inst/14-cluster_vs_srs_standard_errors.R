library(tidyverse)
library(lavaan.bingof)
library(furrr)
plan(multisession, workers = parallel::detectCores() - 1)
theme_set(theme_bw())

# Same population... 
# 1. SRS
# 2. cluster sample
# check power 

n_sim <- 1000
model_no <- 1
samp_size <- 1000
true_vals <- get_true_values(model_no, arrange = c("lambda", "tau", "rho"))
pop <- make_population(model_no, seed = 123)

make_srs_clust_compare <- function(i = 1) {
  fit_srs <- lavaan::sem(
    model = txt_mod(model_no), 
    data = gen_data_bin_srs(pop, n = samp_size), 
    estimator = "PML",
    std.lv = TRUE
  )
  fit_clust <- lavaan::sem(
    model = txt_mod(model_no), 
    data = gen_data_bin_clust(pop, n = samp_size), 
    estimator = "PML",
    std.lv = TRUE, 
    sampling.weights = "wt"
  )      
  se0 <- lavaan::parTable(fit_srs) |> filter(free > 0)  # PML
  se1 <- lavaan::parTable(fit_clust) |> filter(free > 0)  # PML (weighted)
  
  tab <- tibble(
    name  = names(lavaan::coef(fit_srs)),
    kind  = gsub("\\d", "", names(true_vals)),
    truth = true_vals,
    n     = samp_size
  )
  
  res <- 
    bind_cols(
      i = i,
      bind_rows(
        bind_cols(tab, tibble(method = "SRS",
                              est    = se0$est,
                              se     = se0$se)),
        bind_cols(tab, tibble(method = "Clust",
                              est    = se1$est,
                              se     = se1$se))
      )
    )
  
  res
}

res <- 
  future_map(
    .x = 1:n_sim,
    .f = make_srs_clust_compare,
    .options = furrr_options(seed = TRUE),
    .progress = TRUE
  )

p_compare_se_srs_clust <-
  bind_rows(res) |>
  mutate(
    bias = est - truth,
    kind = case_when(
      kind == "lambda" ~ "Loadings",
      kind == "tau"    ~ "Thresholds"
    )
  ) |>
  ggplot(aes(bias, name, fill = method)) +
  geom_boxplot(coef = Inf) +
  geom_vline(linetype = "dashed", xintercept = 0, col = "grey50") +
  facet_grid(kind ~ ., scales = "free") +
  labs(y = NULL, x = "Bias", fill = "Sampling\nmethod") +
  theme_bw()

save(p_compare_se_srs_clust, file = "simulations/p_compare_se_srs_clust.Rdata")


