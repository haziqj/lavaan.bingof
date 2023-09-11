# ORDER OF THE PROBABILITIES IS IMPORTANT. SHOULD JUST USE LAVAAN's ORDERING?
# But possible that not all response patterns are observed...

#' Create a table of binary response patterns
#'
#' @description Ordering is incremental patterns of power 2 from right to left.
#'
#' @param p (integer > 0) The number of items.
#'
#' @return A [tibble()] containing ordinal binary values (0/1) for the items as
#'   well as a column indicating one of the \eqn{2^p} response patterns.
#' @export
#'
#' @examples
#' create_resp_pattern(p = 3)
create_resp_pattern <- function(p = 3) {
  vars <- list()
  for (P in 1:p) vars[[P]] <- 1:0
  names(vars) <- paste0("y", p:1)
  tab <- expand.grid(vars)[p:1] %>%
    as_tibble() %>%
    unite("pattern", everything(), sep = "", remove = FALSE) %>%
    mutate(across(starts_with("y"), ordered)) %>%
    select(starts_with("y"), everything())

  # attr(tab, "patterns") <-
  #   tab %>%
  #   mutate(across(everything(), as.numeric)) %>%
  #   unite("r", everything(), sep = "") %>%
  #   unlist()

  tab
}

#' Create transformation matrices
#'
#' @name transformation-matrices
#' @rdname transformation-matrices
#' @description The derived limited information test statistics involves some
#'   design matrices which act as transformations from the larger \eqn{2^p}
#'   response pattern space to the lower order univariate and bivariate
#'   marginals.
#'
#' - `create_G_mat()` returns the \eqn{\tilde R \times R} indicator matrix to obtain all pairwise components.
#'
#' - `create_T2_mat()` returns the \eqn{p(p+1)/2 \times 2^p} indicator matrix \eqn{T_2} to pick out the unviariate and bivariate moments from the response patterns.
#'
#' - `create_Beta_mat()` returns the \eqn{4p \times p(p+1)/2} design matrix \eqn{\Beta} described in the manuscript (used to express parameters in terms of residuals).
#'
#' Note that ordering is similar to the ordering in [create_resp_pattern()].
#' These design matrices currently only apply to binary data. See technical
#' documents for more details.
#'
#' @inheritParams create_resp_pattern
#'
#' @return A matrix. Additionally, we may inspect the attributes regarding the
#'   ordering of the pairwise components of the \eqn{G} matrix.
#'
#' @examples
#' create_G_mat(p = 3)
#' create_T2_mat(p = 3)
#' create_Beta_mat(p = 3)
NULL

#' @rdname transformation-matrices
#' @export
create_G_mat <- function(p = 3) {
  dat <- create_resp_pattern(p = p) %>%
    select(starts_with("y"))
  id <- combn(p, 2)
  counter <- 0
  G <- list()
  for (k in seq_len(ncol(id))) {
    i <- id[1, k]  # var1
    j <- id[2, k]  # var2
    for (y2 in 0:1) {
      for (y1 in 0:1) {
        # the sequence is 00, 10, 01, 11
        counter <- counter + 1
        G[[counter]] <- as.numeric(dat[, i, drop = TRUE] ==
                                     y1 & dat[, j, drop = TRUE] == y2)
      }
    }
    the_ind <- (counter - 3):counter
    G[[the_ind[1]]] <- - (G[[the_ind[2]]] + G[[the_ind[3]]] + G[[the_ind[4]]])
  }
  G <- as.data.frame(G)
  colnames(G) <- NULL
  G <- as.matrix(G)
  attr(G, "patterns") <- attr(dat, "patterns")
  attr(G, "pairwise") <- id

  t(G)
}

#' @rdname transformation-matrices
#' @export
create_T2_mat <- function(p = 3) {
  dat <- create_resp_pattern(p = p) %>%
    select(starts_with("y"))
  pp <- choose(p, 2)

  M1 <- dat %>%
    mutate(across(everything(), function(x) as.numeric(x) - 1)) %>%
    as.matrix() %>%
    t()

  # Only works for binary data right now... the M2 matrix is picking out every
  # 4th row because for binary data there are total of 4 possible choices within
  # the pairwise groups.
  M2 <- create_G_mat(p = p)[seq(4, 4 * pp, by = 4), ]

  M <- rbind(M1, M2)
  rownames(M) <- NULL
  M
}

#' @rdname transformation-matrices
#' @export
#' @author Myrsini Katsikatsou (`create_Beta_mat()`)
create_Beta_mat <- function(p = 3) Beta_mat_design(nvar = p)

Beta_mat_design <- function(nvar) {
  # Warning: Only applies to binary data.
  # Written by Myrsini Katsikatsou

  no_biv  <- nvar * (nvar - 1) / 2
  idx_univ <- 1:nvar
  idx_biv <- (nvar + 1):(nvar + no_biv)
  idx_cat_var1 <- rep(c(0, 1, 0, 1), times = no_biv)
  idx_cat_var2 <- rep(c(0, 0, 1, 1), times = no_biv)
  idx_pairs <- combn(nvar, 2)
  idx_var1 <- rep(idx_pairs[1, ], each = 4)
  idx_var2 <- rep(idx_pairs[2, ], each = 4)

  nrow_B <- 4 * no_biv
  idx_row <- 1:nrow_B
  B_mat <- matrix(0, nrow = nrow_B, ncol = (nvar + no_biv))

  # Determine the pxixj11, i.e. P(yi = 1, yj = 1)
  diag(B_mat[idx_cat_var1 == 1 & idx_cat_var2 == 1, idx_biv]) <- 1

  # Determine the pxixj10, i.e. P(yi = 1, yj = 0)
  cdtn10 <- idx_cat_var1 == 1 & idx_cat_var2 == 0
  idc_row10 <- idx_row[cdtn10]
  B_mat[cbind(idc_row10, idx_univ[idx_var1[cdtn10]])] <- 1
  diag(B_mat[idc_row10, idx_biv]) <- -1

  # Determine the pxixj01, i.e. P(yi = 0, yj = 1)
  cdtn01 <- idx_cat_var1 == 0 & idx_cat_var2 == 1
  idc_row01 <- idx_row[cdtn01]
  B_mat[cbind(idc_row01, idx_univ[idx_var2[cdtn01]])] <- 1
  diag(B_mat[idc_row01, idx_biv]) <- -1

  # Determine the pxixj00, , i.e. P(yi = 0, yj = 0)
  cdtn00 <- idx_cat_var1 == 0 & idx_cat_var2 == 0
  cdtn11 <- idx_cat_var1 == 1 & idx_cat_var2 == 1
  B_mat[cdtn00, ] <- -1 * (B_mat[cdtn10, ] + B_mat[cdtn01, ] + B_mat[cdtn11, ])

  B_mat
}
