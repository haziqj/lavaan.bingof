# Testing unnormalised weights -------------------------------------------------
dat <- gen_data_bin_wt(1, n = 1000)
tmp <- get_Hinv_mat(select(dat, -wt), 1, dat$wt)


fit <- lavaan::sem(model = txt_mod(1), data = dat, estimator = "PML",
                   std.lv = TRUE, sampling.weights = "wt",
                   sampling.weights.normalization = "none")
# > coef(fit)
# eta1=~y1 eta1=~y2 eta1=~y3 eta1=~y4 eta1=~y5    y1|t1    y2|t1    y3|t1    y4|t1    y5|t1
# 0.792    0.662    0.499    0.332    0.362   -1.415   -0.514   -0.113   -0.790   -1.246
all_tests(fit)


fit <- lavaan::sem(model = txt_mod(1), data = dat, estimator = "PML",
                   std.lv = TRUE, sampling.weights = "wt")
tmp <- calc_test_stuff(fit)
