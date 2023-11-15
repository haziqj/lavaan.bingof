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

# Using manual PL function
fit_man <- optim(theta0, lavaan.bingof:::pl_fn, method = "BFGS",
                 control = list(trace = 10),
                 hessian = TRUE, model_no = model_no, data = dat, wt = NULL)

# Compare coeficients
coef(fit_lav)
fit_man$par

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

get_J_mat <- function(.theta, lavobject = fit_lav) {
  # Code obtained from lavaan:::lav_model_h1_information_firstorder() found in
  # lav_model_h1_information.R
  SIGMA <- lavaan.bingof:::theta_to_Vy(.theta)
  MU    <- rep(0, nrow(SIGMA))
  TH    <- .theta[grepl("tau", names(.theta))]
  PI    <- NULL
  EXO   <- NULL

  lavmodel <- lavobject@Model
  lavcache <- lavobject@Cache
  lavdata <- lavobject@Data
  if(.hasSlot(lavdata, "weights")) {
    WT <- lavdata@weights[[1]]
  } else {
    WT <- NULL
  }

  # Fix lavmodel!!!
  lavmodel@GLIST <- fix_GLIST(.theta, lavobject)
  Delta <- lavaan:::computeDelta(lavmodel)[[1]]

  SC <- lavaan:::pml_deriv1(
    Sigma.hat  = SIGMA,
    Mu.hat     = MU,
    TH         = TH,
    th.idx     = lavmodel@th.idx[[1]],
    num.idx    = lavmodel@num.idx[[1]],
    X          = lavdata@X[[1]],
    eXo        = EXO,
    wt         = NULL,
    PI         = PI,
    lavcache   = lavcache[[1]],
    missing    = lavdata@missing,
    scores     = TRUE,
    negative   = FALSE
  )

  if (is.null(WT)) {
    B1 <- lavaan:::lav_matrix_crossprod(SC)
  } else {
    B1 <- crossprod(WT * SC)
  }

  t(Delta) %*% B1 %*% Delta / 1000
}


HHH <- function(theta) {
  H <- hessian(lavaan.bingof:::pl_fn, theta, model_no = model_no, data = dat)
  1000 * MASS::ginv(H)
}

AAA <- function(theta) {
  tmp <- function(x) {
    H <- HHH(x)  # this is Hinv
    J <- get_J_mat(x)
    -0.5 * sum(diag(H %*% J))
  }
  numDeriv::grad(tmp, theta)
}

theta_hat <- fit_man$par
A <- AAA(theta_hat)
Hinv <- lavaan.bingof:::get_sensitivity_inv_mat(fit_lav)
theta_tilde <- theta_hat + Hinv %*% A
theta_tilde
