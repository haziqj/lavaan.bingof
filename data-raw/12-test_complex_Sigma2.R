library(tidyverse)
library(lavaan.bingof)
library(survey)
library(srvyr)

dat <- gen_data_bin_strcl(population = make_population(5), n = 1000)
svy <- svydesign(ids = ~ school + class, strata = ~ type,
                 weights = ~ wt, data = dat, nest = TRUE)
fit <- lavaan::sem(model = txt_mod(5), data = dat, estimator = "PML",
                   std.lv = TRUE, sampling.weights = "wt")
# Sigma2 <- lavaan.bingof:::create_Sigma2_matrix(fit)
Sigma2 <- lavaan.bingof:::create_Sigma2_matrix_complex(fit, svy)

p <- 15
y_form <- paste0("~ ", paste0("y", 1:p, collapse = " + "), sep = "")
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


my_fun <- function(.lavobject = fit, .svy_design = svy) {
  list2env(extract_lavaan_info(.lavobject), environment())
  list2env(get_uni_bi_moments(.lavobject), environment())

  y_form <- paste0("~ ", paste0("y", 1:p, collapse = " + "), sep = "")
  v <-
    srvyr::as_survey(.svy_design) %>%
    mutate(across(starts_with("y"), \(y) as.numeric(y) - 1))
  ystart <- min(grep("^y[0-9]", names(v$variables))) - 1
  idx <- combn(p, 2)
  for (k in seq_len(ncol(idx))) {
    i <- idx[1, k] + ystart
    j <- idx[2, k] + ystart
    varname <- paste0("y", i, j, collapse = "")
    yi <- v$variables[, i, drop = TRUE]
    yj <- v$variables[, j, drop = TRUE]
    yij <- (yi == 1) * (yj == 1)  # both positive
    v$variables[[varname]] <- yij
    y_form <- paste0(y_form, paste0(" + ", varname))
  }

  x <- v$variables[, -(1:ystart)]
  xbar <- pi2 <- c(pidot1, pidot2)

  x <- t(t(x) - xbar)
  S <- ncol(x)

  a <- matrix(rep(x, S), ncol = S ^ 2)
  b <- x[, rep(1:S, each = S)]

  pweights <- 1 / .svy_design$prob
  psum <- sum(pweights)  # equals N

  cov.wt(x, pweights, center = FALSE)$cov

}























