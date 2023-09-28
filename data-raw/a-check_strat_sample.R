library(tidyverse)
library(lavaan.bingof)
model_no <- 1
run_ligof_sims(1, samp_size = 10000, nsim = 250, samp = "clust")
