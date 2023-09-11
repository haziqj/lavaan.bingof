# population should have 50 strata of unequal sizes Nh. Total pop size is roughly N = 1,000,000.
# strat sample from each strata, select say 20 individuals. The prob of selection is 20 / Nh
#
# make_population2 <- function(model_no = 1, seed = 123, H1 = FALSE,
#                              return_all = FALSE) {
#
#   set.seed(seed)
#
#   # Population settings --------------------------------------------------------
#   # create the population per
#   Nh <- abs(rnorm(50, mean = 1e6 / 50, sd = 1e6 / 50))
#   Nh <- round(1e6 * Nh / sum(Nh), 0)
#   Nh[length(Nh)] <- 1e6 - sum(Nh[-length(Nh)])
#
#   # Set up the loadings and covariance matrices --------------------------------
#   Lambda      <- loading_mat(model_no)
#   neta        <- ncol(Lambda)  # q
#   nitems      <- nrow(Lambda)  # p
#   Psi         <- cov_lv_mat(model_no)
#   Theta       <- matrix(0, nrow = nitems, ncol = nitems)
#   diag(Theta) <- 1 - diag(Lambda %*% Psi %*% t(Lambda))
#   tau         <- get_tau(model_no)
#
#   # Prepare population data ----------------------------------------------------
#   letters <- c(letters, sapply(letters, function(x) paste0(x, letters)))
#   LETTERS <- c(LETTERS, sapply(LETTERS, function(x) paste0(x, LETTERS)))
#
#   pop <-
#     tibble(strata = factor(LETTERS[1:50], levels = LETTERS),
#            Nh     = Nh) %>%
#     rowwise() %>%
#     mutate(type = list(rep(.data$strata, Nh))) %>%
#     unnest_longer(type) %>%
#     select(type)
#   N <- nrow(pop)
#
#   # Generate the data ----------------------------------------------------------
#   eta     <- mvnfast::rmvn(n = N, mu = rep(0, neta), sigma = Psi)
#   epsilon <- mvnfast::rmvn(n = N, mu = rep(0, nitems), sigma = Theta)
#   # eta     <- mnormt::rmnorm(n = N, mean = rep(0, neta), varcov = Psi)
#   # epsilon <- mnormt::rmnorm(n = N, mean = rep(0, nitems), varcov = Theta)
#   ystar   <- tcrossprod(eta, Lambda) + epsilon
#
#   if (isTRUE(H1)) {
#     # Add an extra factor to misspecify the model fit (for power simulations)
#     extra_Lambda <- Lambda[, 1, drop = FALSE] + rnorm(nitems, sd = 0.1)
#     if (model_no <= 3) { extra_Lambda[seq(2, nitems, by = 2), 1] <- 0 }
#     ystar <- ystar + t(extra_Lambda %*% rnorm(N))
#     ystar <- scale(ystar)
#   }
#
#   # repair eta
#   eta <- as.data.frame(eta)
#   colnames(eta) <- paste0("eta", seq_len(neta))
#
#   # Stratify according to latent variable --------------------------------------
#   abil_order <-
#     as_tibble(eta) %>%
#     mutate(z = rowSums(across(everything())),
#            rn = row_number()) %>%
#     select(z, rn) %>%
#     arrange(desc(z)) %>%
#     mutate(type = pop$type) %>%
#     group_by(type) %>%
#     mutate(rn2 = sample(dplyr::n(), replace = FALSE)) %>%
#     arrange(type, rn2) %>%
#     pull(rn)
#   ystar <- ystar[abil_order, ]
#
#   # Get the response patterns --------------------------------------------------
#   y <-
#     t(apply(ystar, 1, function(x) as.numeric(x > tau))) %>%
#     as.data.frame()
#   colnames(y) <- paste0("y", seq_len(nitems))
#   Sigma2 <- (N - 1) / N * cov(convert_dat_to_unibiv(y))
#
#   ystar <- as.data.frame(ystar)
#   colnames(ystar) <- paste0("ystar", seq_len(nitems))
#
#   if (isTRUE(return_all)) {
#     res <- bind_cols(pop, y, ystar, eta[abil_order, , drop = FALSE])
#   } else {
#     res <- bind_cols(pop, y)
#   }
#
#   attr(res, "Sigma2") <- Sigma2
#   res
# }

gen_data_bin_strat2 <- function(model_no = 1, seed = NULL, H1 = FALSE,
                                return_all = FALSE, n = 1000) {
  set.seed(123)
  N <- 1e6
  Nstr <- 50
  Nh <- abs(rnorm(Nstr, mean = N / Nstr, sd = N / Nstr))
  Nh <- round(N * Nh / sum(Nh), 0)
  Nh[length(Nh)] <- N - sum(Nh[-length(Nh)])
  npsu <- floor(n / Nstr)

  seeds <- Nh + seed
  if (is.null(seed)) {
    dats <- purrr::map(Nh, \(x) gen_data_bin(model_no, n = npsu, seed = NULL))
  } else {
    dats <- purrr::map2(Nh, seeds,
                        \(x, y) gen_data_bin(model_no, n = npsu, seed = y))
  }

  tibble(
    stratum = seq_len(Nstr),
    wt = Nh / npsu,
    dat = dats
  ) %>%
    tidyr::unnest_wider(dat) %>%
    tidyr::unnest_longer(starts_with("y")) %>%
    mutate(wt = wt / sum(wt) * dplyr::n())
}
