library(lavaan)
# library(lavaan.bingof)
library(tidyverse)

M2_sim <- function(i) {
  # Code to generate 5 item CFA data
  popmod <- "
  f1 =~ 0.8*y1 + 0.7*y2 + 0.5*y3 + 0.4*y4 + 0.3*y5
  f1 ~~ 1*f1

  y1| -1.43*t1
  y2| -0.55*t2
  y3| -0.13*t3
  y4| -0.72*t4
  y5| -1.13*t5
"
  dat <- simulateData(popmod, sample.nobs = 2500)
  dat <- as.data.frame(lapply(dat, ordered))

  mod <- "f1 =~ y1 + y2 + y3 + y4 + y5"
  fit <- cfa(mod, dat, std.lv = TRUE, estimator = "PML")

  # Ingredients for the ligof test
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

  n <- nobs(fit)
  Delta2 <- create_Delta2_matrix(fit)
  Delta2comp <- mcompanion::null_complement(Delta2)
  Sigma2 <- lavaan.bingof:::create_Sigma2_matrix(fit)
  with(lavaan.bingof:::get_uni_bi_moments(fit), {
    p2_hat <<- c(pdot1, pdot2)
    pi2_hat <<- c(pidot1, pidot2)
  })
  e2_hat <- p2_hat - pi2_hat

  C2 <-
    Delta2comp %*%
    MASS::ginv(t(Delta2comp) %*% Sigma2 %*% Delta2comp) %*%
    t(Delta2comp)

  M2 <- n * colSums(e2_hat * (C2 %*% e2_hat))
  df <- nrow(Delta2) - ncol(Delta2)
  pval <- pchisq(M2, df, lower.tail = FALSE)

  tibble(M2 = M2, df = df, pval = pval)
}

library(furrr)
plan(multisession)

res <- future_map(1:10000, M2_sim, .progress = TRUE, .options = furrr_options(seed = TRUE))

bind_rows(res) |>
  ggplot(aes(sample = M2)) +
  stat_qq(distribution = qchisq, dparams = list(df = 5)) +
  stat_qq_line(distribution = qchisq, dparams = list(df = 5)) +
  labs(title = "QQ-plot of M2 against chi-square distribution with 5 df") +
  theme_bw()

bind_rows(res) |>
  ggplot(aes(x = M2, y = ..density..)) +
  geom_histogram(col = "gray90") +
  geom_line(
    data = tibble(x = seq(0, 20, 0.1), y = dchisq(x, df = 5)),
    aes(x, y),
    col = "red3",
    linewidth = 1
  ) +
  theme_bw()

bind_rows(res) |>
  ggplot(aes(x = pval, y = ..density..)) +
  geom_histogram(binwidth = 0.05, boundary = 0, col = "white") +
  geom_hline(yintercept = 1, col = "red3", linetype = "dashed") +
  theme_bw()


