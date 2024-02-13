library(tidyverse)
library(lavaan.bingof)
library(lavaan)
library(numDeriv)

model_no <- 1
mod <- txt_mod(model_no)
dat <- gen_data_bin(model_no)
theta0 <- get_true_values(model_no)

# Using lavaan
fit_lav <- sem(mod, dat, std.lv = TRUE, estimator = "PML")
theta_hat <- theta0
theta_hat[] <- coef(fit_lav)

# Hinv <- lavaan.bingof:::get_sensitivity_inv_mat(fit_lav_fxd)
# lavaan:::lav_model_nvcov_robust_sandwich
# lavaan:::lav_model_information_firstorder(fit_lav_fxd)
# lavaan:::lav_model_information_firstorder(
#   lavmodel = fit_lav_fxd@Model,
#   lavsamplestats = fit_lav_fxd@SampleStats,
#   lavdata = fit_lav_fxd@Data,
#   lavcache = fit_lav_fxd@Cache,
#   lavimplied = fit_lav_fxd@implied,
#   lavh1 = fit_lav_fxd@h1,
#   lavoptions = fit_lav_fxd@Options,
#   extra = FALSE,
#   check.pd = FALSE,
#   augmented = FALSE,
#   inverted = TRUE,
#   use.ginv = TRUE
# )

# Fix param values
mod_fix <- function(.theta) {
  # For model 1 only!!!
  res <-"
  eta1 =~ lambda1*y1 + lambda2*y2 + lambda3*y3 + lambda4*y4 + lambda5*y5

  y1| tau1*t1
  y2| tau2*t1
  y3| tau3*t1
  y4| tau4*t1
  y5| tau5*t1
  "
  for (i in seq_along(.theta)) {
    res <- gsub(names(.theta)[i], .theta[i], res)
  }
  res
}

fix_GLIST <- function(theta = theta0, lavobject = fit_lav) {
  lambdas <- theta[grepl("lambda", names(theta))]
  rhos <- theta[grepl("rho", names(theta))]
  taus <- theta[grepl("tau", names(theta))]

  GLIST <- lavobject@Model@GLIST
  GLIST$lambda[] <- lambdas
  GLIST$theta <- diag(1 - diag(tcrossprod(GLIST$lambda)))
  GLIST$tau[] <- taus
  GLIST
}

# Function to grab the J matrix
JJJ <- function(theta = theta_hat, lavobject = fit_lav) {
  # Fit a fixed parameter model
  fit_lav_fxd <- sem(mod_fix(theta), dat, std.lv = TRUE, estimator = "PML")
  lavmodel <- lavobject@Model
  lavmodel@GLIST <- fix_GLIST(theta, lavobject)
  Delta <- lavaan:::computeDelta(lavmodel)[[1]]
  B1 <- lavaan:::lav_model_h1_information_firstorder(
    lavmodel = fit_lav_fxd@Model,
    lavsamplestats = fit_lav_fxd@SampleStats,
    lavdata = fit_lav_fxd@Data,
    lavoptions = fit_lav_fxd@Options,
    lavimplied = fit_lav_fxd@implied,
    lavh1 = fit_lav_fxd@h1,
    lavcache = fit_lav_fxd@Cache
  )[[1]]

  # Return the unit information matrix
  t(Delta) %*% B1 %*% Delta / lavobject@Data@nobs[[1]]
}

# get_J_mat <- function(.theta, lavobject = fit_lav) {
#   # Code obtained from lavaan:::lav_model_h1_information_firstorder() found in
#   # lav_model_h1_information.R
#   SIGMA <- lavaan.bingof:::theta_to_Vy(.theta)
#   MU    <- rep(0, nrow(SIGMA))
#   TH    <- .theta[grepl("tau", names(.theta))]
#   PI    <- NULL
#   EXO   <- NULL
#
#   lavmodel <- lavobject@Model
#   lavcache <- lavobject@Cache
#   lavdata <- lavobject@Data
#   if(.hasSlot(lavdata, "weights")) {
#     WT <- lavdata@weights[[1]]
#   } else {
#     WT <- NULL
#   }
#
#   # Fix lavmodel!!!
#   lavmodel@GLIST <- fix_GLIST(.theta, lavobject)
#   Delta <- lavaan:::computeDelta(lavmodel)[[1]]
#
#   SC <- lavaan:::pml_deriv1(
#     Sigma.hat  = SIGMA,
#     Mu.hat     = MU,
#     TH         = TH,
#     th.idx     = lavmodel@th.idx[[1]],
#     num.idx    = lavmodel@num.idx[[1]],
#     X          = lavdata@X[[1]],
#     eXo        = EXO,
#     wt         = NULL,
#     PI         = PI,
#     lavcache   = lavcache[[1]],
#     missing    = lavdata@missing,
#     scores     = TRUE,
#     negative   = FALSE
#   )
#
#   if (is.null(WT)) {
#     B1 <- lavaan:::lav_matrix_crossprod(SC)
#   } else {
#     B1 <- crossprod(WT * SC)
#   }
#
#   t(Delta) %*% B1 %*% Delta / 1000
# }

# Using manual PL function
fit_man <- optim(theta0, lavaan.bingof:::pl_fn, method = "BFGS",
                 control = list(trace = 10),
                 hessian = TRUE, model_no = model_no, data = dat, wt = NULL)

# Compare coeficients
coef(fit_lav)
fit_man$par

HHH <- function(theta = theta_hat) {
  H <- hessian(lavaan.bingof:::pl_fn, theta, model_no = model_no, data = dat)
  H / nrow(dat)
}

AAA <- function(.theta) {
  tmp <- function(x) {
    Hinv <- MASS::ginv(HHH(x))
    J    <- JJJ(x)
    -0.5 * sum(diag(Hinv %*% J))
  }
  numDeriv::grad(tmp, .theta)
}

A <- AAA(theta_hat)
Hinv <- lavaan.bingof:::get_sensitivity_inv_mat(fit_lav)
theta_tilde <- theta_hat + Hinv %*% A
theta_tilde
