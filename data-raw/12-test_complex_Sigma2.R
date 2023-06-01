library(tidyverse)
library(lavaan.bingof)
library(survey)
library(srvyr)

dat <- gen_data_bin_strcl(n = 1000)
svy <- svydesign(ids = ~ school + class, strata = ~ type,
                 weights = ~ wt, data = dat, nest = TRUE)
fit <- lavaan::sem(model = txt_mod(1), data = dat, estimator = "PML",
                   std.lv = TRUE, sampling.weights = "wt")
Sigma2 <- lavaan.bingof:::create_Sigma2_matrix(fit)

p <- 5
y_form <- paste0("~ ", paste0("y", 1:5, collapse = " + "), sep = "")
v <-
  as_survey(svy) %>%
  mutate(across(starts_with("y"), as.numeric))
ystart <- min(grep("^y[0-9]", names(v$variables))) - 1
idx <- combn(p, 2) + ystart
for (k in seq_len(ncol(idx))) {
  i <- idx[1, k]
  j <- idx[2, k]
  varname <- paste0("y", i - ystart, j - ystart, collapse = "")
  yi <- v$variables[, i, drop = TRUE]
  yj <- v$variables[, j, drop = TRUE]
  yij <- 1 + (yi == 2) * (yj == 2)  # both positive
  v$variables[[varname]] <- yij
  y_form <- paste0(y_form, paste0(" + ", varname))
}
v <- v %>%
  svyvar(as.formula(y_form), .) %>%
  as.matrix()
attr(v, "var") <- NULL
attr(v, "statistic") <- NULL
