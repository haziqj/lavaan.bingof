library(tidyverse)
library(lavaan.bingof)
library(lme4)
library(furrr)
plan(multisession, workers = parallel::detectCores() - 1)

# For all models, check bias in the strcl sampling
repeat_strcl_sampling <- function(model_no = 1, samp_size = 1000, B = 1000) {
  
  cli::cli_alert_info("Model {model_no}, sample size {samp_size}...")
  pop <- lavaan.bingof:::make_population(model_no)
  
  my_fun <- function(i) {
    dat <- gen_data_bin_strcl(pop, n = samp_size)
    fit <- lavaan::sem(
      model = txt_mod(model_no), 
      data = dat, 
      estimator = "PML",
      std.lv = TRUE, 
      sampling.weights = "wt"
    )
    true_vals <- get_true_values(model_no, arrange = c("lambda", "tau", "rho"))
    se <- lavaan::parTable(fit) |> filter(free > 0)
    tibble(
      i     = i,
      name  = names(lavaan::coef(fit)),
      kind  = gsub("\\d", "", names(true_vals)),
      truth = true_vals,
      n     = samp_size,
      est   = se$est,
      se    = se$se
    )
  }
  
  furrr::future_map(
    .x = 1:B,
    .f = possibly(my_fun, NA),
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
  )
}

res <-
  expand_grid(
    model_no = 1:5,
    samp_size = c(500, 1000, 2500, 5000, 10000)
  ) |>
  mutate(
    res = map2(model_no, samp_size, repeat_strcl_sampling)
  )

tab_strcl_bias <-
  res |>
  mutate(res = map(res, bind_rows)) |>
  unnest(res) |>
  summarise(
    bias = mean(est - truth),
    .by = c(samp_size, model_no)
  ) |>
  pivot_wider(
    id_cols = samp_size,
    names_from = model_no,
    values_from = bias
  )

save(tab_strcl_bias, file = "simulations/tab_strcl_bias.RData")
