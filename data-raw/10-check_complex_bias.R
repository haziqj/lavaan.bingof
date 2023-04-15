library(tidyverse)
library(ggridges)
library(kableExtra)
library(lavaan)
library(survey)
library(truncnorm)
theme_set(theme_bw())
library(doSNOW)
library(foreach)
library(lavaan.bingof)

pop <- make_population(1)
true_vals <- get_true_values(1)
mod <- txt_mod(1)

B <- 100  # no of simulations
ncores <- parallel::detectCores() - 2
pb <- txtProgressBar(min = 0, max = B, style = 3)
progress <- function(i) setTxtProgressBar(pb, i)
cl <- makeCluster(ncores)
registerDoSNOW(cl)

res <- foreach(
  b = 1:B,
  .combine = bind_rows,
  # .errorhandling = "remove",
  .packages = c("lavaan", "tidyverse", "lavaan.bingof"),
  .options.snow = list(progress = progress)
) %dopar% {

  # Stratified sampling
  dat <- gen_data_bin_complex1(population = pop)
  fit11 <- sem(model = mod, data = dat, estimator = "PML", std.lv = TRUE)
  fit21 <- sem(model = mod, data = dat, estimator = "PML", std.lv = TRUE,
               sampling.weights = "wt")

  # 2-stage cluster sampling
  dat <- gen_data_bin_complex2(population = pop)
  fit12 <- sem(model = mod, data = dat, estimator = "PML", std.lv = TRUE)
  fit22 <- sem(model = mod, data = dat, estimator = "PML", std.lv = TRUE,
               sampling.weights = "wt")

  # 2-stage stratified cluster sampling
  dat <- gen_data_bin_complex3(population = pop)
  fit13 <- sem(model = mod, data = dat, estimator = "PML", std.lv = TRUE)
  fit23 <- sem(model = mod, data = dat, estimator = "PML", std.lv = TRUE,
               sampling.weights = "wt")

  bind_rows(
    tibble(
      B = b, param = names(coef(fit11)), Sampling = "Stratified",
      type = c(rep("Loadings", 5), rep("Thresholds", 5)),
      "truth" = as.numeric(true_vals),
      "PML" = as.numeric(coef(fit11)), "PML (wt)" = as.numeric(coef(fit21))
    ),
    tibble(
      B = b, param = names(coef(fit11)), Sampling = "2S Cluster",
      type = c(rep("Loadings", 5), rep("Thresholds", 5)),
      "truth" = as.numeric(true_vals),
      "PML" = as.numeric(coef(fit12)), "PML (wt)" = as.numeric(coef(fit22))
    ),
    tibble(
      B = b, param = names(coef(fit11)), Sampling = "2S Strat-Clust",
      type = c(rep("Loadings", 5), rep("Thresholds", 5)),
      "truth" = as.numeric(true_vals),
      "PML" = as.numeric(coef(fit13)), "PML (wt)" = as.numeric(coef(fit23))
    )
  )
}
close(pb)
stopCluster(cl)

save(res, file = "check_complex_bias.RData")
