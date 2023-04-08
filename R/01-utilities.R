# This script contains the utility functions for data generation, including
# lavaan textual models, getting the true parameter values, etc.

## ---- Libraries --------------------------------------------------------------
library(tidyverse)
# remotes::install_github("haziqj/lavaan")  # modified lavaan to do PL weights
library(lavaan)
library(mnormt)
library(doSNOW)
library(foreach)

## ---- Textual model ----------------------------------------------------------
txt_mod <- function(model.no) {
  # Function to return the textual model for lavaan fit, depending on which
  # simulation scenario we are investigating.
  if (model.no == 1) mod <- "eta1 =~ NA*y1 + y2 + y3 + y4 + y5"
  if (model.no == 2) mod <- "eta1 =~ NA*y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8"
  if (model.no == 3) mod <- "eta1 =~ NA*y1 +  y2 +  y3 +  y4 +  y5 +
                                        y6 +  y7 +  y8 +  y9 + y10 +
                                       y11 + y12 + y13 + y14 + y15"
  if (model.no == 4) mod <- "eta1 =~ NA*y1 + y2 + y3 + y4 + y5
                             eta2 =~ NA*y6 + y7 + y8 + y9 + y10"
  if (model.no == 5) mod <- "eta1 =~  NA*y1 +  y2 +  y3 +  y4 +  y5
                             eta2 =~  NA*y6 +  y7 +  y8 +  y9 + y10
                             eta3 =~ NA*y11 + y12 + y13 + y14 + y15"

  return(mod)
}

## ---- Loading matrix ---------------------------------------------------------
loading_mat <- function(model.no) {
  lam_entries <- c(0.80, 0.70, 0.47, 0.38, 0.34)
  lam_index <- rep(seq_along(lam_entries), 3)

  if (model.no == 1) { nitems <- 5;  neta <- 1 }
  if (model.no == 2) { nitems <- 8;  neta <- 1 }
  if (model.no == 3) { nitems <- 15; neta <- 1 }
  if (model.no == 4) { nitems <- 10; neta <- 2 }
  if (model.no == 5) { nitems <- 15; neta <- 3 }

  res <- matrix(0, nrow = nitems, ncol = neta)

  if (model.no %in% 1:3) {
    res[, 1] <- lam_entries[lam_index[seq_len(nitems)]]
  } else {
    for (k in seq_len(neta)) {
      res[(k - 1) * 5 + 1:5, k] <- lam_entries
    }
  }

  return(res)
}

get_Lambda <- loading_mat  # alternative name

## ---- Cov LV matrix ----------------------------------------------------------
cov_lv_mat <- function(model.no) {
  if (model.no %in% 1:3) {
    return(1)
  }
  if (model.no == 4) {
    return(matrix(c(1, 0.3, 0.3, 1), nrow = 2, ncol = 2))
  }
  if (model.no == 5) {
    return(matrix(c(  1, 0.2, 0.3,
                    0.2,   1, 0.4,
                    0.3, 0.4, 1), nrow = 3, ncol = 3))
  }
}

get_Psi <- cov_lv_mat

## ---- Get thresholds----------------------------------------------------------
get_tau <- function(model.no = 1) {
  nitems <- nrow(loading_mat(model.no))
  tau <- rep(c(-1.43, -0.55, -0.13, -0.72, -1.13), 3)
  tau[seq_len(nitems)]
}

## ---- Get corr matrix --------------------------------------------------------
get_Sigmay <- function(model.no) {
  Lambda      <- loading_mat(model.no)
  Phi         <- cov_lv_mat(model.no)
  neta        <- ncol(Lambda)
  nitems      <- nrow(Lambda)
  Theta       <- matrix(0, nrow = nitems, ncol = nitems)
  diag(Theta) <- 1 - diag(Lambda %*% Phi %*% t(Lambda))

  Lambda %*% Phi %*% t(Lambda) + Theta
}

## ---- Get true values --------------------------------------------------------
get_true_values <- function(model.no) {
  # Gets the true values of the freely estimated theta values in this order:
  # lambda (loadings), rho (factor correlations), tau (thresholds)

  # Loadings -------------------------------------------------------------------
  Lambda <- get_Lambda(model.no)
  p <- nrow(Lambda)
  q <- ncol(Lambda)
  lambda <- c(Lambda)
  # if (q > 1) {
  #   names(lambda) <- paste0(
  #     "lambda",
  #     matrix(1:p, nrow = p, ncol = q),
  #     ",",
  #     matrix(1:q, nrow = p, ncol = q, byrow = TRUE)
  #   )
  # } else {
  #   names(lambda) <- paste0("lambda", seq_len(p))
  # }
  lambda <- lambda[lambda != 0]
  names(lambda) <- paste0("lambda", seq_along(lambda))

  # Factor correlations --------------------------------------------------------
  Psi <- get_Psi(model.no)
  rho <- c(Psi[upper.tri(Psi)])
  q <- nrow(Psi)
  if (!is.null(q)) {
    names(rho) <- expand.grid(1:q, 1:q)[c(upper.tri(Psi)), ] %>%
      apply(., 1, paste0, collapse = "") %>%
      paste0("rho", .)
  }

  # Threholds ------------------------------------------------------------------
  tau <- get_tau(model.no)
  names(tau) <- paste0("tau", seq_along(tau))

  c(lambda, rho, tau)

}
