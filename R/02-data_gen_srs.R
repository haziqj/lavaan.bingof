#' Generate independent data samples according to `model_no`
#'
#' @inheritParams txt_mod
#' @param n (numeric > 0) Sample size.
#' @param seed (numeric) The random seed.
#' @param H1 (logical) Are we generating a data to misspecify the true model?
#'   For power simulations.
#' @param return_all (logical) Return the underlying latent variable \eqn{y^*} as
#'   well?
#'
#' @return A [tibble()] containing ordinal binary values (0/1) for the items.
#' @export
#'
#' @examples
#' gen_data_bin(1)
gen_data_bin <- function(model_no = 1, n = 1000, seed = NULL, H1 = FALSE,
                         return_all = FALSE) {
  # Generate data for simple random sample
  set.seed(seed)
  orig.model_no <- model_no

  # Set up the loadings and covariance matrices --------------------------------
  Lambda      <- loading_mat(model_no)
  neta        <- ncol(Lambda)
  nitems      <- nrow(Lambda)
  Psi         <- cov_lv_mat(model_no)
  Theta       <- matrix(0, nrow = nitems, ncol = nitems)
  diag(Theta) <- 1 - diag(Lambda %*% Psi %*% t(Lambda))
  tau         <- get_tau(model_no)

  # Generate the data ----------------------------------------------------------
  eta   <- mnormt::rmnorm(n = n, mean = rep(0, neta), varcov = Psi)
  delta <- mnormt::rmnorm(n = n, mean = rep(0, nitems), varcov = Theta)
  ystar <- t(Lambda %*% t(eta)) + delta

  if (isTRUE(H1)) {
    # Add an extra factor to misspecify the model fit (for power simulations)
    extra_Lambda <- Lambda[, 1, drop = FALSE] + rnorm(nitems, sd = 0.1)
    if (model_no <= 3) { extra_Lambda[seq(2, nitems, by = 2), 1] <- 0 }
    ystar <- ystar + t(extra_Lambda %*% rnorm(n))
    ystar <- scale(ystar)
  }

  y <-
    t(apply(ystar, 1, function(x) as.numeric(x > tau))) %>%
    as.data.frame() %>%
    mutate(across(everything(), ordered))
  colnames(y) <- paste0("y", seq_len(nitems))

  ystar <- as.data.frame(ystar)
  colnames(ystar) <- paste0("ystar", seq_len(nitems))

  if (isTRUE(return_all)) {
    bind_cols(y, ystar)
  } else {
    as_tibble(y)
  }
}
