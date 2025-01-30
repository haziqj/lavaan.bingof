library(lavaan)
library(lavaan.bingof)
library(tidyverse)
library(tinytest)

popmod <- "
  f1 =~ 0.8*y1 + 0.7*y2 + 0.5*y3 + 0.4*y4 + 0.3*y5
  f1 ~~ 1*f1

  y1| -1.43*t1
  y2| -0.55*t2
  y3| -0.13*t3
  y4| -0.72*t4
  y5| -1.13*t5
"
dat <- simulateData(popmod, sample.nobs = 1000)
dat <- as.data.frame(lapply(dat, ordered))

mod <- "f1 =~ y1 + y2 + y3 + y4 + y5"
fit <- cfa(mod, dat, std.lv = TRUE, estimator = "PML")

## ----- New Sigma2 function ---------------------------------------------------

create_Delta2_matrix <- function(lavobject) {

  p <- lavobject@Model@nvar
  all_thresh <- inspect(lavobject, "est")$tau

  if(ncol(all_thresh) != 1L) {
    stop("This simplified function only handles purely binary indicators (1 threshold per variable).")
  }
  tau <- as.numeric(all_thresh)
  Sigma_hat <- inspect(lavobject, "implied")$cov
  rho_ij <- Sigma_hat[lower.tri(Sigma_hat)]

  Delta_full <- lavaan:::computeDelta(lavobject@Model)[[1]]
  derTauToTheta <- Delta_full[1:p, , drop = FALSE]  # dtau/dtheta
  derRhoToTheta <- Delta_full[-(1:p), , drop = FALSE]  #drho_ij/dtheta

  # Mini function for bivariate derivatives
  get_dP11 <- function(taui, tauj, rho) {
    dP.dTaui <- -1 * dnorm(taui) * pnorm((rho * taui - tauj) / sqrt(1 - rho^2))
    dP.dTauj <- -1 * dnorm(tauj) * pnorm((rho * tauj - taui) / sqrt(1 - rho^2))
    dP.dRho  <- mvtnorm::dmvnorm(c(taui, tauj), sigma = matrix(c(1, rho, rho, 1), 2, 2))

    c(dP.dTaui = dP.dTaui, dP.dTauj = dP.dTauj, dP.dRho  = dP.dRho)
  }

  # Now compute bivariate derivatives
  pair_idx <- which(lower.tri(Sigma_hat), arr.ind = TRUE)
  npairs <- nrow(pair_idx)
  derBiv_11_wrtTheta <- matrix(0, nrow = npairs, ncol = ncol(derTauToTheta))

  for (idx in seq_len(npairs)) {
    i <- pair_idx[idx, 2]
    j <- pair_idx[idx, 1]
    dP_local <- get_dP11(tau[i], tau[j], Sigma_hat[i, j])

    # Chain rule for tau_i, tau_j
    dP_taui_theta <- dP_local["dP.dTaui"] * derTauToTheta[i, ]
    dP_tauj_theta <- dP_local["dP.dTauj"] * derTauToTheta[j, ]
    dP_rho_theta <- dP_local["dP.dRho"] * derRhoToTheta[idx, ]

    derBiv_11_wrtTheta[idx, ] <- dP_taui_theta + dP_tauj_theta + dP_rho_theta
  }

  # Return
  rbind(
    -1 * dnorm(tau) * derTauToTheta,  # derUni_1_wrtTheta
    derBiv_11_wrtTheta
  )

}

create_Delta2_matrix_optimised <- function(lavobject) {
  p <- lavobject@Model@nvar
  all_thresh <- inspect(lavobject, "est")$tau

  if(ncol(all_thresh) != 1L) {
    stop("This simplified function only handles purely binary indicators (1 threshold per variable).")
  }
  tau <- as.numeric(all_thresh)
  Sigma_hat <- inspect(lavobject, "implied")$cov
  rho_ij <- Sigma_hat[lower.tri(Sigma_hat)]

  Delta_full <- lavaan:::computeDelta(lavobject@Model)[[1]]
  derTauToTheta <- Delta_full[1:p, , drop = FALSE]
  derRhoToTheta <- Delta_full[-(1:p), , drop = FALSE]

  pair_idx <- which(lower.tri(Sigma_hat), arr.ind = TRUE)
  npairs <- nrow(pair_idx)

  # Precompute all necessary values
  dnorm_tau <- dnorm(tau)
  tau_i <- tau[pair_idx[, 2]]
  tau_j <- tau[pair_idx[, 1]]
  dnorm_tau_i <- dnorm_tau[pair_idx[, 2]]
  dnorm_tau_j <- dnorm_tau[pair_idx[, 1]]
  rho_values <- rho_ij

  # Common terms for vectorized calculations
  denominator_sq <- 1 - rho_values^2
  sqrt_denominator_sq <- sqrt(denominator_sq)
  z1 <- (rho_values * tau_i - tau_j) / sqrt_denominator_sq
  z2 <- (rho_values * tau_j - tau_i) / sqrt_denominator_sq

  # Vectorized derivatives calculations
  dP.dTaui <- -dnorm_tau_i * pnorm(z1)
  dP.dTauj <- -dnorm_tau_j * pnorm(z2)
  exponent <- -0.5 * (tau_i^2 - 2 * rho_values * tau_i * tau_j + tau_j^2) / denominator_sq
  dP.dRho <- exp(exponent) / (2 * pi * sqrt_denominator_sq)

  # Vectorized matrix operations
  i_indices <- pair_idx[, 2]
  j_indices <- pair_idx[, 1]

  dP_taui_theta <- dP.dTaui * derTauToTheta[i_indices, , drop = FALSE]
  dP_tauj_theta <- dP.dTauj * derTauToTheta[j_indices, , drop = FALSE]
  dP_rho_theta <- dP.dRho * derRhoToTheta

  derBiv_11_wrtTheta <- dP_taui_theta + dP_tauj_theta + dP_rho_theta

  # Combine results
  rbind(
    -dnorm_tau * derTauToTheta,
    derBiv_11_wrtTheta
  )
}

Delta2 <- create_Delta2_matrix(fit)
Delta2op <- create_Delta2_matrix_optimised(fit)
Delta2lb <- lavaan.bingof:::get_Delta_mats(fit)$Delta2
tinytest::expect_equal(Delta2, Delta2lb)
tinytest::expect_equal(Delta2, Delta2op)
tinytest::expect_equal(Delta2op, Delta2lb)

