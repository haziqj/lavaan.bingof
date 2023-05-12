library(tidyverse)
library(lavaan.bingof)
library(foreach)
library(doSNOW)

# Initialise parallel stuff
nsim <- 100
ncores <- 30
pb <- txtProgressBar(min = 0, max = nsim, style = 3)
progress <- function(i) setTxtProgressBar(pb, i)
cl <- makeCluster(ncores)
registerDoSNOW(cl)

res <- foreach(
  i = 1:nsim, #.combine = bind_rows,
  .packages = c("lavaan.bingof", "tidyverse", "lavaan", "survey"),
  .export = ls(globalenv()),
  .errorhandling = "remove",
  .options.snow = list(progress = progress)
) %dopar% {
  fit <- lavaan.bingof:::fit_facmod_pml(1, "strat", n = 3000, H1 = TRUE)
  all_tests(fit$fit, fit$svy)
}

bind_rows(res) %>%
  mutate(rej = pval < 0.05) %>%
  group_by(name) %>%
  summarise(rate = mean(rej))

dat <- gen_data_bin_strat(n = 3000)
