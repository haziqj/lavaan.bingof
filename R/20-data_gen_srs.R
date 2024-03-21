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
  neta        <- ncol(Lambda)  # q
  nitems      <- nrow(Lambda)  # p
  Psi         <- cov_lv_mat(model_no)
  Theta       <- matrix(0, nrow = nitems, ncol = nitems)
  diag(Theta) <- 1 - diag(Lambda %*% Psi %*% t(Lambda))
  tau         <- get_tau(model_no)

  # Generate the data ----------------------------------------------------------
  eta     <- mvnfast::rmvn(n = n, mu = rep(0, neta), sigma = Psi)
  epsilon <- mvnfast::rmvn(n = n, mu = rep(0, nitems), sigma = Theta)
  # eta   <- mnormt::rmnorm(n = n, mean = rep(0, neta), varcov = Psi)
  # delta <- mnormt::rmnorm(n = n, mean = rep(0, nitems), varcov = Theta)
  extra <- 0

  if (isTRUE(H1)) {
    # Add an extra factor to misspecify the model fit (for power simulations)
    extra <- rnorm(n)
    Lambda <- cbind(Lambda, apply(Lambda, 1, sum))
    if (model_no %in% c(1, 2)) {
      rem <- c(4, 5)
      Lambda[rem, 1] <- Lambda[-rem, 2] <- 0
    }
    if (model_no == 3) {
      rem <- c(3, 4, 5)
      Lambda[rem, 1] <- Lambda[-rem, 2] <- 0
    }
    if (model_no == 4) {
      Lambda[1, 1] <- 0
      Lambda[nitems, 2] <- 0
      Lambda[-c(1, nitems), 3] <- 0
    }
    if (model_no == 5) {
      rem <- c(3, 8)
      Lambda[rem[1], 1] <- 0
      Lambda[rem[2], 2] <- 0
      # Lambda[rem[3], 3] <- 0
      Lambda[-rem, 4] <- 0
    }
  }

  ystar <- tcrossprod(cbind(eta, extra), Lambda) + epsilon
  y <-
    {1 * (ystar > matrix(tau, nrow = n, ncol = nitems, byrow = TRUE))} |>
    as.data.frame() |>
    mutate(across(everything(), \(x) ordered(x, levels = c(0, 1))))
  colnames(y) <- paste0("y", seq_len(nitems))

  ystar <- as.data.frame(ystar)
  colnames(ystar) <- paste0("ystar", seq_len(nitems))

  if (isTRUE(return_all)) {
    eta <- as.data.frame(eta)
    colnames(eta) <- paste0("eta", seq_len(neta))
    as_tibble(bind_cols(y, ystar, eta))
  } else {
    as_tibble(y)
  }
}
