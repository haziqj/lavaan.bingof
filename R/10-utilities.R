#' Function to return the textual model for lavaan fit
#'
#' @param model_no (integer) Choose from 1--5. See pkgdown articles for details.
#'
#' @return Character vector of factor model to input to [lavaan::lavaan()].
#' @export
#'
#' @examples
#' txt_mod(1)
#' txt_mod(2)
#' txt_mod(3)
#' txt_mod(4)
#' txt_mod(5)
txt_mod <- function(model_no = 1L) {

  if (model_no == 1) mod <- "eta1 =~ NA*y1 + y2 + y3 + y4 + y5"
  if (model_no == 2) mod <- "eta1 =~ NA*y1 + y2 + y3 + y4 + y5 + y6 + y7 + y8"
  if (model_no == 3) mod <- "eta1 =~ NA*y1 +  y2 +  y3 +  y4 +  y5 +
                                        y6 +  y7 +  y8 +  y9 + y10 +
                                       y11 + y12 + y13 + y14 + y15"
  if (model_no == 4) mod <- "eta1 =~ NA*y1 + y2 + y3 + y4 + y5
                             eta2 =~ NA*y6 + y7 + y8 + y9 + y10"
  if (model_no == 5) mod <- "eta1 =~  NA*y1 +  y2 +  y3 +  y4 +  y5
                             eta2 =~  NA*y6 +  y7 +  y8 +  y9 + y10
                             eta3 =~ NA*y11 + y12 + y13 + y14 + y15"

  mod
}

## ---- Loading matrix ---------------------------------------------------------
loading_mat <- function(model_no) {
  lam_entries <- rep(c(0.80, 0.70, 0.47, 0.38, 0.34), 3)
  # lam_entries <- c(0.80, 0.70, 0.60, 0.50, 0.40,
  #                  0.85, 0.75, 0.65, 0.55, 0.42,
  #                  0.90, 0.77, 0.67, 0.57, 0.45)

  if (model_no == 1) { nitems <- 5;  neta <- 1 }
  if (model_no == 2) { nitems <- 8;  neta <- 1 }
  if (model_no == 3) { nitems <- 15; neta <- 1 }
  if (model_no == 4) { nitems <- 10; neta <- 2 }
  if (model_no == 5) { nitems <- 15; neta <- 3 }

  res <- matrix(0, nrow = nitems, ncol = neta)

  if (model_no %in% 1:3) {
    res[, 1] <- lam_entries[seq_len(nitems)]
  } else {
    for (k in seq_len(neta)) {
      res[(k - 1) * 5 + 1:5, k] <- lam_entries[(k - 1) * 5 + 1:5]
    }
  }

  res
}

get_Lambda <- loading_mat  # alternative name

## ---- Cov LV matrix ----------------------------------------------------------
cov_lv_mat <- function(model_no) {
  if (model_no %in% 1:3) {
    return(1)
  }
  if (model_no == 4) {
    return(matrix(c(1, 0.3, 0.3, 1), nrow = 2, ncol = 2))
  }
  if (model_no == 5) {
    return(matrix(c(  1, 0.2, 0.3,
                    0.2,   1, 0.4,
                    0.3, 0.4, 1), nrow = 3, ncol = 3))
  }
}

get_Psi <- cov_lv_mat

## ---- Get thresholds----------------------------------------------------------
get_tau <- function(model_no = 1) {
  nitems <- nrow(loading_mat(model_no))
  tau <- rep(c(-1.43, -0.55, -0.13, -0.72, -1.13), 3)
  # tau <- c(-1.5, -0.9, -0.3, 0.3, 0.9,
  #          -1.3, -0.7, -0.1, 0.5, 1.1,
  #          -1.1, -0.5, 0.1, 0.7, 1.3)
  tau[seq_len(nitems)]
}

## ---- Get corr matrix --------------------------------------------------------
get_Sigmay <- function(model_no) {
  Lambda      <- loading_mat(model_no)
  Phi         <- cov_lv_mat(model_no)
  neta        <- ncol(Lambda)
  nitems      <- nrow(Lambda)
  Theta       <- matrix(0, nrow = nitems, ncol = nitems)
  diag(Theta) <- 1 - diag(Lambda %*% Phi %*% t(Lambda))

  Lambda %*% Phi %*% t(Lambda) + Theta
}

#' Returns theoretical true values
#'
#' @description
#'
#' - `get_true_values()` returns the true values of the freely estimated theta values in this order: \eqn{\lambda} (loadings), \eqn{\rho} (factor correlations), \eqn{\tau} (thresholds).
#'
#' - `get_theoretical_uni_bi_moments()` returns the univariate (`pidot1`) and bivariate (`pidot2`) theoretical probabilities of successes.
#'
#' @inheritParams txt_mod
#' @param arrange How should the true values be arranged? By default it is in the order of loadings, factor correlations, and thresholds.
#'
#' @return A vector of true parameter values used for the simulations.
#' @export
#'
#' @seealso [gen_data_bin()], [gen_data_bin_strat()],
#'   [gen_data_bin_clust()], [gen_data_bin_strcl()], [gen_data_bin_srs()]
#'
#' @examples
#' get_true_values(1)
get_true_values <- function(model_no, arrange = c("lambda", "rho", "tau")) {

  # Loadings -------------------------------------------------------------------
  Lambda <- get_Lambda(model_no)
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
  Psi <- get_Psi(model_no)
  rho <- c(Psi[upper.tri(Psi)])
  q <- nrow(Psi)
  if (!is.null(q)) {
    names(rho) <- expand.grid(1:q, 1:q)[c(upper.tri(Psi)), ] %>%
      apply(., 1, paste0, collapse = "") %>%
      paste0("rho", .)
  }

  # Thresholds -----------------------------------------------------------------
  tau <- get_tau(model_no)
  names(tau) <- paste0("tau", seq_along(tau))

  res <- list(lambda = lambda, rho = rho, tau = tau)
  res <- res[arrange]
  names(res) <- NULL
  unlist(res)

}

#' @rdname get_true_values
#' @param collapse (logical) Should a vector be returned instead of a list
#'   separating the univariate and bivariate quantities?
#' @inherit get_uni_bi_moments return
#' @export
#' @examples
#' get_theoretical_uni_bi_moments(1)
get_theoretical_uni_bi_moments <- function(model_no, collapse = FALSE) {
  Var_ystar <- get_Sigmay(model_no)
  mu_ystar <- rep(0, nrow(Var_ystar))
  TH <- get_tau(model_no)
  p <- nrow(get_Lambda(model_no))

  # Univariate -----------------------------------------------------------------
  pidot1 <- rep(NA, p)
  for (i in seq_along(pidot1)) {
    pidot1[i] <- pnorm(TH[i], mean = mu_ystar[i], sd = sqrt(Var_ystar[i, i]),
                       lower.tail = FALSE)
  }

  # Bivariate ------------------------------------------------------------------
  id <- combn(p, 2)
  pidot2 <- rep(NA, ncol(id))
  for (k in seq_along(pidot2)) {
    i <- id[1, k]  # var1
    j <- id[2, k]  # var2
    pidot2[k] <- mnormt::sadmvn(lower = c(TH[i], TH[j]), upper = c(Inf, Inf),
                                mean = mu_ystar[c(i, j)],
                                varcov = Var_ystar[c(i, j), c(i, j)])
  }

  if (isTRUE(collapse)) {
    c(pidot1, pidot2)
  } else {
    list(pidot1 = pidot1, pidot2 = pidot2)
  }
}

#' Get univariate and bivariate moments
#'
#' @description Returns univariate and bivariate moments (i.e. positive
#'   probabilities only) based on model i.e. `pidot1` and `pidot2` and
#'   (weighted) sample i.e. `pdot1` and `pdot2`.
#'
#'
#' @param .lavobject A [lavaan::lavaan()] fit object.
#' @param wtd (logical) Should the weighted proportions be used?
#'
#' @returns A list of univariate and bivariate moments.
#' @export
#'
#' @seealso [get_theoretical_uni_bi_moments()]
#'
#' @examples
#' fit <- lavaan::sem(txt_mod(1), gen_data_bin(1, n = 500), std.lv = TRUE,
#'                    estimator = "PML")
#' get_uni_bi_moments(fit)
get_uni_bi_moments <- function(.lavobject, wtd = TRUE) {
  list2env(extract_lavaan_info(.lavobject), environment())
  if (!isTRUE(wtd)) wt <- 1
  N <- sum(wt)

  # Univariate -----------------------------------------------------------------
  pdot1 <- pidot1 <- rep(NA, p)
  for (i in seq_along(pidot1)) {
    pdot1[i] <- sum(wt[dat[, i] == 2]) / N
    pidot1[i] <- pnorm(TH[i], mean = mu_ystar[i], sd = sqrt(Var_ystar[i, i]),
                       lower.tail = FALSE)
  }

  # Bivariate ------------------------------------------------------------------
  id <- combn(p, 2)
  pdot2 <- pidot2 <- rep(NA, ncol(id))
  for (k in seq_along(pidot2)) {
    i <- id[1, k]  # var1
    j <- id[2, k]  # var2
    pdot2[k] <- sum(wt[dat[, i] == 2 & dat[, j] == 2]) / N
    pidot2[k] <- mnormt::sadmvn(lower = c(TH[i], TH[j]), upper = c(Inf, Inf),
                                mean = mu_ystar[c(i, j)],
                                varcov = Var_ystar[c(i, j), c(i, j)])
  }

  list(
    # Univariate moments
    pdot1  = pdot1,   # sample
    pidot1 = pidot1,  # model

    # Bivariate moments
    pdot2 = pdot2,    # sample
    pidot2 = pidot2   # model
  )
}

# Useful function to help fit any of the models we want to analyse -------------
fit_facmod_pml <- function(model_no, samp = c("srs", "wtd", "strat", "clust",
                                              "strcl", "strat2"),
                           n = 1000, seed = NULL, H1 = FALSE,
                           ignore_weights = FALSE) {
  # Convenience function to fit one of our 5 models using lavaan's PML estimator
  # (with weights if necessary) just by specifying the model number. Mostly used
  # for testing so will suppress warning messages.
  samp <- match.arg(samp, c("srs", "wtd", "strat", "clust", "strcl", "strat2"))
  seed_used <- seed
  the_wt <- NULL
  if (samp == "srs") {
    # Simple random sampling ---------------------------------------------------
    dat <- gen_data_bin(model_no = model_no, n = n, seed = seed_used, H1 = H1)
  } else if (samp == "wtd") {
    # Informative sampling -----------------------------------------------------
    dat <- gen_data_bin_wt(model_no = model_no, n = n, seed = seed_used,
                           H1 = H1)
  } else {
    the_wt <- "wt"
    if (samp == "strat2") {
      # Stratified sampling v2 -------------------------------------------------
      dat <- gen_data_bin_strat2(model_no, n = n, seed = seed_used, H1 = H1)
    } else {
      pop <- make_population(model_no, seed = seed, H1 = H1)
      if (samp == "strat") {
        # Stratified sampling --------------------------------------------------
        dat <- gen_data_bin_strat(population = pop, n = n, seed = seed_used)
      }
      if (samp == "clust") {
        # Cluster sampling -----------------------------------------------------
        dat <- gen_data_bin_clust(population = pop, n = n, seed = seed_used)
      }
      if (samp == "strcl") {
        # Stratified-cluster sampling ------------------------------------------
        dat <- gen_data_bin_strcl(population = pop, n = n, seed = seed_used)
      }
    }
  }
  if (isTRUE(ignore_weights)) the_wt <- NULL

  suppressWarnings(
    fit <- lavaan::sem(model = txt_mod(model_no), data = dat, estimator = "PML",
                       std.lv = TRUE, sampling.weights = the_wt, do.fit = TRUE)
  )

  list(dat = dat, fit = fit)
}

neff <- function(x) {
  ntilde <- sum(x) ^ 2 / sum(x ^ 2)
  cli::cli_alert("Effective sample size: {round(ntilde, 2)} / {length(x)} ({round(100 * ntilde/length(x), 2)}%)")
}
