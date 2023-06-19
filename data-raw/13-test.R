res <- run_ligof_sims(3, nsim = 10, samp = "strcl")

library(survey)
library(svrep)

model_no <- 5
dat <- gen_data_bin_strcl(make_population(model_no))
svy <- svydesign(ids = ~ school + class, strata = ~ type, weights = ~ wt,
                 data = dat, nest = TRUE)
fit <- lavaan::sem(model = txt_mod(model_no), data = dat, estimator = "PML",
                   std.lv = TRUE, sampling.weights = "wt")
# bootsvy <- as_bootstrap_design(svy, type = "Rao-Wu-Yue-Beaumont",
#                                replicates = 1000)
all_tests(fit, svy, bootstrap = TRUE)




res <- lavaan.bingof:::fit_facmod_pml(1, "strcl")
res <- lavaan.bingof:::fit_facmod_pml(3, "strcl")
