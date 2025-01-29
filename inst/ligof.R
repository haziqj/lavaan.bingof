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
  dat <- simulateData(popmod, sample.nobs = 1000)
  dat <- as.data.frame(lapply(dat, ordered))

  mod <- "f1 =~ y1 + y2 + y3 + y4 + y5"
  fit <- cfa(mod, dat, std.lv = TRUE, estimator = "PML")

  # Ingredients for the ligof test
  create_Delta2_matrix <- function(lavobject) {

    p <- lavobject@Model@nvar
    all_thresh <- inspect(fit, "est")$tau

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
      dP_local <- get_dP11(tau[i], tau[j], Rho_hat[i, j])

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

res <- future_map(1:1000, M2_sim, .progress = TRUE, .options = furrr_options(seed = TRUE))

# qqplot of M2 against chisquare distribution with 5 df
bind_rows(res) |>
  ggplot(aes(sample = M2)) +
  stat_qq(distribution = qchisq, dparams = list(df = 5)) +
  stat_qq_line(distribution = qchisq, dparams = list(df = 5)) +
  labs(title = "QQ-plot of M2 against chi-square distribution with 5 df") +
  theme_bw()



