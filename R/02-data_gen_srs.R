gen_data_bin <- function(model.no = 1, n = 1000, seed = NULL, H1 = FALSE,
                         return_all = FALSE) {
  # Generate data for simple random sample
  set.seed(seed)

  # Set up the loadings and covariance matrices --------------------------------
  Lambda      <- loading_mat(model.no)
  neta        <- ncol(Lambda)
  nitems      <- nrow(Lambda)
  Psi         <- cov_lv_mat(model.no)
  Theta       <- matrix(0, nrow = nitems, ncol = nitems)
  diag(Theta) <- 1 - diag(Lambda %*% Psi %*% t(Lambda))
  tau         <- get_tau(model.no)

  # Generate the data ----------------------------------------------------------
  eta   <- mnormt::rmnorm(n = n, mean = rep(0, neta), varcov = Psi)
  delta <- mnormt::rmnorm(n = n, mean = rep(0, nitems), varcov = Theta)
  ystar <- t(Lambda %*% t(eta)) + delta

  if (isTRUE(H1)) {
    # Add an extra factor to misspecify the model fit (for power simulations)
    ystar <- ystar +
      t((Lambda[, 1, drop = FALSE] + rnorm(nitems, sd = 0.1)) %*% rnorm(n))
    ystar <- scale(ystar)
  }

  y <-
    t(apply(ystar, 1, function(x) as.numeric(x > tau))) %>%
    as.data.frame() %>%
    mutate_all(ordered)
  colnames(y) <- paste0("y", seq_len(nitems))

  ystar <- as.data.frame(ystar)
  colnames(ystar) <- paste0("ystar", seq_len(nitems))

  if (isTRUE(return_all)) {
    bind_cols(y, ystar)
  } else {
    as_tibble(y)
  }
}
