res <- run_ligof_sims(1, nsim = 250, samp = "strcl")

library(survey)
library(svrep)

dat <- gen_data_bin_strcl(make_population(1))
svy <- svydesign(ids = ~ school + class, strata = ~ type, weights = ~ wt,
                 data = dat, nest = TRUE)
bootsvy <- as_bootstrap_design(svy, type = "Rao-Wu-Yue-Beaumont",
                               replicates = 1000)
