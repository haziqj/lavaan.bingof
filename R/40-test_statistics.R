#' @importFrom MASS ginv
#' @importFrom Matrix rankMatrix
#' @importFrom gtools combinations
#' @importFrom mnormt sadmvn

# Export lavaan functions
lav_model_vcov <- utils::getFromNamespace("lav_model_vcov", "lavaan")
LongVecTH.Rho <- utils::getFromNamespace("LongVecTH.Rho", "lavaan")
computeDelta <- utils::getFromNamespace("computeDelta", "lavaan")
dbinorm <- utils::getFromNamespace("dbinorm", "lavaan")
lav_tables_oneway <- utils::getFromNamespace("lav_tables_oneway", "lavaan")
lav_tables_pairwise_model_pi <- utils::getFromNamespace("lav_tables_pairwise_model_pi", "lavaan")
lav_tables_pairwise_freq_cell <- utils::getFromNamespace("lav_tables_pairwise_freq_cell", "lavaan")
univariateExpProbVec <- utils::getFromNamespace("univariateExpProbVec", "lavaan")

# No visible binding notes
globalVariables(c(".",".data", "Delta2", "N", "Omega2", "S", "Sigma2", "TH",
                  "X2", "Var_ystar", "avg_class_size", "dat", "e2_hat",
                  "lavcache", "lavdata", "lavmodel", "lavoptions",
                  "lavpartable", "lavsamplestats", "mu_ystar", "nschools",
                  "nstudents", "p", "pdot1", "pdot2", "pi2_hat", "pidot1",
                  "pidot2", "pr_class_selected", "pr_school_selected", "prob",
                  "rn", "rn2", "school", "students_in_school_type", "th.idx",
                  "type", "var1", "var2", "wt", "z", "Omega2_rank", "Sigmahat",
                  "alpha_", "converged", "fit", "mean_X2", "mean_df", "name",
                  "pval", "rej_rate", "approx_delta"))

## ---- Utilities --------------------------------------------------------------
extract_lavaan_info <- function(lavobject) {
  lavdata        <- lavobject@Data
  lavmodel       <- lavobject@Model
  lavsamplestats <- lavobject@SampleStats
  lavoptions     <- lavobject@Options
  lavcache       <- lavobject@Cache
  lavpartable    <- lavobject@ParTable
  # lavtable       <- suppressWarnings(lavaan::lavTables(lavobject, dimension = 0))

  TH      <- lavobject@Fit@TH[[1]]
  th.idx  <- lavobject@Model@th.idx[[1]]
  if(lavobject@Model@nexo == 0) {
    nvar <- lavobject@Model@nvar
  } else {
    stop("Case with covariates has not been done yet.")
  }

  # Distribution of ystar (underlying variables)
  Sigmahat <- lavobject@Fit@Sigma.hat[[1]]  #g for multigroup analysis
  Meanhat <- c(lavobject@implied$mean[[1]]) #g for multigroup analysis
  nsize <- lavobject@Data@nobs[[1]]

  # Raw data tibble
  dat <- as.data.frame(lavdata@X) %>% as_tibble()
  colnames(dat) <- lavdata@ov.names[[1]]
  wt <- lavdata@weights[[1]]
  if (is.null(wt)) wt <- rep(1, nsize)
  dat$wt <- wt

  list(
    lavdata        = lavdata,  # data saved in lavaan object
    lavmodel       = lavmodel,
    lavsamplestats = lavsamplestats,
    lavoptions     = lavoptions,
    lavcache       = lavcache,
    lavpartable    = lavpartable,
    TH             = TH,  # thresholds
    th.idx         = th.idx,   # threshold index
    p              = nvar,  # no. of items
    N              = nsize,
    dat            = dat,
    Var_ystar      = Sigmahat,  # varcov matrix of underlying MVN dist
    mu_ystar       = Meanhat  # mean vector of underlying MVN dist
  )
}

get_sensitivity_inv_mat <- function(.lavobject, matrix_type = c("Sensitivity",
                                                                "Sandwich")) {
  # Returns the matrix H^{-1}, the inverse of the sensitivity matrix
  list2env(extract_lavaan_info(.lavobject), environment())
  matrix_type <- match.arg(matrix_type, c("Sensitivity", "Sandwich"))

  # This extracts the inverse of the observed sensitivity matrix H = E(nabla^2
  # pl). It is stored as E.inv as part of the Godambe information matrix (E.inv
  # %*% B0 %*% E.inv) because lavoptions has se = "robust.huber.white"
  VCOV <- lav_model_vcov(lavmodel       = lavmodel,
                         lavsamplestats = lavsamplestats,
                         lavoptions     = lavoptions,
                         lavdata        = lavdata,
                         lavpartable    = lavpartable,
                         lavcache       = lavcache)
  # Note: VCOV here is (H %*% J^{-1} %*% H)^{-1} = H^{-1} %*% J %*% H^{-1} and
  # it is the "full" information matrix. E.inv is therefore the unit information
  # matrix, and J has the full information.

  if (matrix_type == "Sensitivity") {
    mat <- attr(VCOV, "E.inv")  # not yet divided by sqrt n
  } else if (matrix_type == "Sandwich") {
    mat <- VCOV  # already divided by sqrt n
  }

  tt <- sum(lavpartable$free != 0)
  if (is.null(mat)) mat <- diag(tt)
  mat
}

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

## ---- Delta matrices ---------------------------------------------------------
.derModelUnivBivProbToTheta <- function(nvar, TH, th.idx, Sigmahat, lavcache,
                                        lavmodel) {
  # This function returns a list for the derivatives of the probabilities for
  # use in the Delta matrix. Written by Myrsini Katsikatsou.

  # Set up indices -------------------------------------------------------------
  idxTH_var1 <- lavcache$LONG$index.thres.var1.of.pair
  idxTH_var2 <- lavcache$LONG$index.thres.var2.of.pair
  idxLastTH_var1 <- lavcache$LONG$last.thres.var1.of.pair
  idxLastTH_var2 <- lavcache$LONG$last.thres.var2.of.pair

  idx_pairs_ext <- lavcache$LONG$index.pairs.extended
  cdtn          <- (idxTH_var1 != 0) & (idxTH_var2 != 0)
  idx_pairs_S   <- idx_pairs_ext[cdtn]
  idxVar1_S     <- lavcache$LONG$index.var1.of.pair[cdtn]
  idxVar2_S     <- lavcache$LONG$index.var2.of.pair[cdtn]
  idxCat_var1   <- idxTH_var1[cdtn]
  idxCat_var2   <- idxTH_var2[cdtn]

  cors <- Sigmahat[lower.tri(Sigmahat)]
  th.rho.vec <- LongVecTH.Rho(no.x = nvar,
                              all.thres = TH,
                              index.var.of.thres = th.idx,
                              rho.xixj = cors)

  # Delta below gives the derivatives of standardised thresholds and polychoric
  # correlations (rows of Delta in this order) with respect to model parameter
  # vector theta. In lavaan in a factor analysis model the order of the
  # individual parameters is: loadings, unstandardized thresholds, factor
  # correlations
  Delta <- computeDelta(lavmodel = lavmodel)[[1]]
  noTH <- length(TH)
  derRhoToTheta <- Delta[-c(1:noTH), ]  # der. rhoij wrt theta
  derTauToTheta <- Delta[c(1:noTH), ]  # der. tau (unstd.) wrt theta

  # Section 1: der.pi.xixj.ab.to.Rho %*% der.Rho.to.Theta ----------------------
  prob.vec <- rep(NA, length(idxTH_var1))
  prob.vec[idxTH_var1 == 0 | idxTH_var2 == 0 |
             idxLastTH_var1 | idxLastTH_var2] <- 0
  prob.vec[is.na(prob.vec)] <- dbinorm(th.rho.vec$thres.var1.of.pair,
                                       th.rho.vec$thres.var2.of.pair,
                                       rho = th.rho.vec$rho.vector)

  den.term1 <- prob.vec[idxTH_var1 != 0 & idxTH_var2 != 0]
  den.term2 <- prob.vec[idxTH_var1 != 0 & !idxLastTH_var2]
  den.term3 <- prob.vec[idxTH_var2 != 0 & !idxLastTH_var1]
  den.term4 <- prob.vec[!idxLastTH_var1 & !idxLastTH_var2]

  der.pi.xixj.ab.to.rho.xixj <- den.term1 - den.term2 - den.term3 + den.term4

  # Each row of derRhoToTheta corresponds to one rhoij. Thus it needs to be
  # repeated as many times as the cells of the bivariate frequency table of the
  # corresponding pair of variables (xi,xj) so that it can then be mulitplied
  # with der.pi.xixj.to.rho.xixj. That's why we do the following:
  der.pi.xixj.ab.to.theta_part1 <-
    derRhoToTheta[idx_pairs_S, ] * der.pi.xixj.ab.to.rho.xixj

  # der.pi.xixj.ab.to.theta_part1 has as many rows as the bivariate
  # probabilities pi.xixj.ab and as many columns as the dimension of theta. Note
  # that in pi.xixj.ab, index "a" (index of category for var 1) runs the
  # fastest, index "b" (index of the vategory of variable 2 of the pair) is the
  # next fastest, next it is "j" (index of variable 2 of the pair), and last it
  # is "i" (index of variable 1 of the pair)

  # Section 2: der.pi.xixj.ab.to.tau.xi %*% der.Tau.to.Theta -------------------
  # To compute der.pi.xixj.1b.to.tau.xi, note that der.pi.xixj.to.tau.xi include
  # only those pi.xixj.ab where the category "a" corresponds to a threshold that
  # is free to be estimated. In the case of binary data it means the derivatives
  # of pi.xixj.1b. The derivatives of pi.xixj.2b are (-1)*der.pi.xixj.1b.
  xi <- lapply(lavcache$LONG[c("index.thres.var2.of.pair",
                               "last.thres.var2.of.pair")], function(y){
                                 y[!(idxTH_var1==0 | idxLastTH_var1)]
                               })

  cum.prob.vec <- rep(NA, length(xi$index.thres.var2.of.pair))
  cum.prob.vec[xi$index.thres.var2.of.pair == 0] <- 0
  cum.prob.vec[xi$last.thres.var2.of.pair] <- 1
  denom <- sqrt(1 - (th.rho.vec$rho.vector * th.rho.vec$rho.vector))
  cum.prob.vec[is.na(cum.prob.vec)] <- pnorm(
    (th.rho.vec$thres.var2.of.pair -
       th.rho.vec$rho.vector * th.rho.vec$thres.var1.of.pair) / denom
  )
  den.prob.vec <- dnorm(th.rho.vec$thres.var1.for.dnorm.in.der.pi.to.tau.xi)
  der.pi.xixj.1b.to.tau.xi <- den.prob.vec *
    (cum.prob.vec[ xi$index.thres.var2.of.pair != 0] -
       cum.prob.vec[!xi$last.thres.var2.of.pair])

  # Since the derivatives of pi.xixj.2b are (-1)*der.pi.xixj.1b , the complete
  # vector der.pi.xixj.ab.to.tau.xi is
  idc_der.pi.xixj.1b.to.tau.xi <- rep(TRUE, length(idx_pairs_S))
  idc_der.pi.xixj.1b.to.tau.xi[idxCat_var1 == 2] <- FALSE
  der.pi.xixj.ab.to.tau.xi <- rep(NA, length(idx_pairs_S))
  der.pi.xixj.ab.to.tau.xi[idc_der.pi.xixj.1b.to.tau.xi] <-
    der.pi.xixj.1b.to.tau.xi
  der.pi.xixj.ab.to.tau.xi[!idc_der.pi.xixj.1b.to.tau.xi] <-
    (-1) * der.pi.xixj.1b.to.tau.xi

  # Note that for binary data each row of derTauToTheta corresponds to one tau.i
  # where i the vairable index, while der.pi.xixj.ab.to.tau.xi is a vector of
  # length equal to as many as bivariate variables pi.xixj.ab. Thus each row of
  # derTauToTheta needs to be repeated as many times as the cells of the
  # bivariate frequency table of the corresponding pair of variables (xi,xj) (4
  # cells for binary data). Then, it can then be mulitplied with
  # der.pi.xixj.ab.to.tau.xi. That's why we do the following:
  der.pi.xixj.ab.to.theta_part2 <-
    derTauToTheta[idxVar1_S,] * der.pi.xixj.ab.to.tau.xi

  # Section 3: der.pi.xixj.ab.to.tau.xj %*% der.Tau.to.Theta -------------------
  xj <- lapply(lavcache$LONG[c("index.thres.var1.of.pair",
                               "last.thres.var1.of.pair")],
               function(y){ y[!(idxTH_var2==0 | idxLastTH_var2)] })

  cum.prob.vec <- rep(NA, length(xj$index.thres.var1.of.pair) )
  cum.prob.vec[xj$index.thres.var1.of.pair==0] <- 0
  cum.prob.vec[xj$last.thres.var1.of.pair] <- 1
  cum.prob.vec[is.na(cum.prob.vec)] <- pnorm(
    (th.rho.vec$thres.var1.of.pair -
       th.rho.vec$rho.vector* th.rho.vec$thres.var2.of.pair) / denom
  )                                                          # denom from above
  den.prob.vec <- dnorm(th.rho.vec$thres.var2.for.dnorm.in.der.pi.to.tau.xj)
  der.pi.xixj.a1.to.tau.xj <-
    den.prob.vec * (cum.prob.vec[ xj$index.thres.var1.of.pair!=0] -
                      cum.prob.vec[!xj$last.thres.var1.of.pair])

  # Note that der.pi.xixj.to.tau.xj include only those pi.xixj.ab where the
  # category "b" corresponds to a threshold that is free to be estimated. In the
  # case of binary data only the derivative of pi.xixj.a1 is computed. We use
  # the fact that der.pi.xixj.a2 = - der.pi.xixj.a1.to.tau.xj to complete the
  # vector  der.pi.xixj.ab.to.tau.xj.
  idc_der.pi.xixj.a1.to.tau.xj <- rep(TRUE, length(idx_pairs_S))
  idc_der.pi.xixj.a1.to.tau.xj[idxCat_var2 == 2] <- FALSE
  der.pi.xixj.ab.to.tau.xj <- rep(NA, length(idx_pairs_S))
  der.pi.xixj.ab.to.tau.xj[idc_der.pi.xixj.a1.to.tau.xj] <-
    der.pi.xixj.a1.to.tau.xj
  der.pi.xixj.ab.to.tau.xj[!idc_der.pi.xixj.a1.to.tau.xj] <-
    (-1) * der.pi.xixj.a1.to.tau.xj

  # As said, for binary data, each row of derTauToTheta corresponds to one tau.i
  # where i the variable index, while der.pi.xixj.ab.to.tau.xj is a vector of
  # length equal to as many as bivariate variables pi.xixj.ab. Thus each row of
  # derTauToTheta needs to be repeated as many times as the cells of the
  # bivariate frequency table of the corresponding pair of variables (xi,xj) (4
  # cells for binary data). Then, it can then be mulitplied with
  # der.pi.xixj.ab.to.tau.xj. That's why we do the following:
  der.pi.xixj.ab.to.theta_part3 <-
    derTauToTheta[idxVar2_S, ] * der.pi.xixj.ab.to.tau.xj

  # Section 4: der.pi.xixj.ab.to.theta and der.pi.xixj.22.to.theta (positive ---
  # bivariate probs) -----------------------------------------------------------
  der.pi.xixj.ab.to.theta_F <-
    der.pi.xixj.ab.to.theta_part1 + der.pi.xixj.ab.to.theta_part2 +
    der.pi.xixj.ab.to.theta_part3

  der.pi.xixj.22.to.theta_F <-
    der.pi.xixj.ab.to.theta_F[(idxCat_var1 == 2 & idxCat_var2 == 2), ]

  # Section 5: der.pi.xi.1.to.theta --------------------------------------------
  # Note that for our problem we need ONLY the derivatives of the POSITIVE
  # univariate probabilities w.r.t to theta and NOT of all all univariate
  # probabilities. The univariate probabilities are only functions of the
  # standardised thresholds. The command below is correct only in the case of
  # binary data.
  derUnivProb1wrtTau <- (-1) * dnorm(TH)
  derUnivProb1wrtTheta <- derTauToTheta * derUnivProb1wrtTau

  # Return output --------------------------------------------------------------
  list(derBivProbwrtTheta   = der.pi.xixj.ab.to.theta_F,
       derBivProb11wrtTheta = der.pi.xixj.22.to.theta_F,
       derUnivProb1wrtTheta = derUnivProb1wrtTheta,
       Delta = Delta,
       idx_pairs_ext,
       cdtn,
       idx_pairs_S,
       idxVar1_S,
       idxVar2_S,
       idxCat_var1,
       idxCat_var2)
}

get_Delta_mats <- function(.lavobject) {
  list2env(extract_lavaan_info(.lavobject), environment())

  derivatives <- .derModelUnivBivProbToTheta(
    nvar = p, TH = TH, th.idx = th.idx, Sigmahat = Var_ystar,
    lavcache = lavcache[[1]], lavmodel = lavmodel
  )

  list(
    Delta2    = rbind(derivatives$derUnivProb1wrtTheta,
                      derivatives$derBivProb11wrtTheta),
    Delta_til = derivatives$derBivProbwrtTheta
  )
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
  if (!isTRUE(wtd)) dat$wt <- 1

  # Univariate -----------------------------------------------------------------
  pdot1 <- pidot1 <- rep(NA, p)
  for (i in seq_along(pidot1)) {
    pdot1[i] <- sum(dat$wt[dat[, i] == 2]) / N
    pidot1[i] <- pnorm(TH[i], mean = mu_ystar[i], sd = sqrt(Var_ystar[i, i]),
                       lower.tail = FALSE)
  }

  # Bivariate ------------------------------------------------------------------
  id <- combn(p, 2)
  pdot2 <- pidot2 <- rep(NA, ncol(id))
  for (k in seq_along(pidot2)) {
    i <- id[1, k]  # var1
    j <- id[2, k]  # var2
    pdot2[k] <- sum(dat$wt[dat[, i] == 2 & dat[, j] == 2]) / N
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

## ---- Sigma multinomial matrices ---------------------------------------------
create_Sigma2_matrix <- function(.lavobject) {
  list2env(extract_lavaan_info(.lavobject), environment())
  list2env(get_uni_bi_moments(.lavobject), environment())

  S <- p * (p + 1) / 2
  Eysq <- matrix(NA, S, S)

  id <- c(1:p, asplit(combn(p, 2), 2))
  # idS <- t(combn(S, 2))
  idS <- gtools::combinations(S, 2, repeats = TRUE)
  colnames(idS) <- c("i", "j")
  idy <- as_tibble(idS) %>%
    mutate(var1 = id[i], var2 = id[j],
           y = mapply(c, var1, var2, SIMPLIFY = FALSE))

  for (s in seq_len(nrow(idS))) {
    i <- idy$i[s]
    j <- idy$j[s]
    yy <- unique(idy$y[[s]])

    dimy <- length(yy)
    Eysq[i, j] <- mnormt::sadmvn(lower = TH[yy], upper = rep(Inf, dimy),
                                 mean = mu_ystar[yy], varcov = Var_ystar[yy, yy])

  }
  Eysq[lower.tri(Eysq)] <- t(Eysq)[lower.tri(Eysq)]
  Eysq - tcrossprod(c(pidot1, pidot2))
}

# create_Sigma2_matrix_complex_old <- function(.lavobject, .svy_design) {
#   list2env(extract_lavaan_info(.lavobject), environment())
#   list2env(get_uni_bi_moments(.lavobject), environment())
#
#   # Get strata information
#   strata  <- .svy_design$strata[, 1, drop = FALSE]
#   colnames(strata) <- "strata"
#   a <- unique(strata[, 1])
#
#   # Get cluster information
#   cluster <- .svy_design$cluster
#   colnames(cluster)[1] <- "cluster"
#
#   dat <-
#     dat %>%
#     bind_cols(strata, cluster)
#
#   # Create the data set of univariate and bivariate responses (1/0)
#   tmp <- dat %>%
#     select(starts_with("y")) %>%
#     mutate(across(starts_with("y"), function(x) as.numeric(x == 2)))
#
#   id <- combn(p, 2)
#   the_names <- paste0("y", 1:p)
#   the_names <- c(the_names, paste0(the_names[id[1, ]], the_names[id[2, ]]))
#   tmp2 <- as_tibble(((tmp[, id[1, ]] == 1) + (tmp[, id[2, ]] == 1)) == 2) %>%
#     mutate(across(everything(), as.numeric)) %>% suppressWarnings()
#   colnames(tmp2) <- the_names[-(1:p)]
#   tmp <- bind_cols(dat %>% select(-starts_with("y")), tmp, tmp2)
#
#   pi2 <- c(pidot1, pidot2)
#   Sigma <- list()
#   for (i in a) {
#     # Filter the strata a
#     current_dat <- tmp %>%
#       filter(strata == i)
#
#     # Obtain the responses in this strata (i.e. according to S_ab)
#     y <- current_dat %>%
#       select(starts_with("y")) %>%
#       mutate(across(starts_with("y"), function(x) as.numeric(x == 1))) %>%
#       as.matrix()
#
#     uab <- {sweep(y, 2, pi2, "-") * current_dat$wt} %>%
#       as_tibble() %>%
#       bind_cols(current_dat %>% select(-starts_with("y")), .) %>%
#       group_by(cluster) %>%
#       summarise(across(starts_with("y"), sum), .groups = "drop") %>%
#       select(starts_with("y")) %>%
#       as.matrix()
#     na <- nrow(uab)
#
#     ubar <- apply(uab, 2, mean)
#
#     Sigma[[i]] <- crossprod(sweep(uab, 2, ubar, "-")) * na / (na - 1)
#   }
#   Reduce("+", Sigma) / N
# }

create_Sigma2_matrix_complex <- function(.lavobject, .svy_design,
                                         bootstrap = FALSE, nboot = 100) {
  list2env(extract_lavaan_info(.lavobject), environment())
  list2env(get_uni_bi_moments(.lavobject), environment())

  y_form <- paste0("~ ", paste0("y", 1:p, collapse = " + "), sep = "")
  v <-
    srvyr::as_survey(.svy_design) %>%
    mutate(across(starts_with("y"), \(y) as.numeric(y) - 1))
  ystart <- min(grep("^y[0-9]", names(v$variables))) - 1
  idx <- combn(p, 2)
  for (k in seq_len(ncol(idx))) {
    i <- idx[1, k]
    j <- idx[2, k]
    varname <- paste0("y", i, ".", j, collapse = "")
    yi <- v$variables[, i + ystart, drop = TRUE]
    yj <- v$variables[, j + ystart, drop = TRUE]
    yij <- (yi == 1) * (yj == 1)  # both positive
    v$variables[[varname]] <- yij
    y_form <- paste0(y_form, paste0(" + ", varname))
  }

  x <- v$variables[, -(1:ystart)]
  xbar <- c(pidot1, pidot2)  # pi2 (model probs)
  # xbar <- c(pdot1, pdot2)  #p2 (proportions)
  # x <- t(t(x) - xbar)

  if (isTRUE(bootstrap)) {
    bootsvy <- svrep::as_bootstrap_design(.svy_design,
                                          type = "Rao-Wu-Yue-Beaumont",
                                          replicates = nboot)
    wt <- bootsvy$repweights
    res <- NULL
    for (b in seq_len(ncol(wt))) {
      res[[b]] <- cov.wt(x, wt[, b], center = FALSE)$cov
    }
    res <- apply(simplify2array(res), 1:2, mean)
  } else {
    wt <- .svy_design$allprob$wt
    res <- cov.wt(x, wt, center = xbar)$cov
  }

  res
}

# create_Sigma2_matrix_complex2 <- function(.lavobject, .svy_design) {
#   list2env(extract_lavaan_info(.lavobject), environment())
#   list2env(get_uni_bi_moments(.lavobject), environment())
#
#   y_form <- paste0("~ ", paste0("y", 1:p, collapse = " + "), sep = "")
#   v <-
#     srvyr::as_survey(.svy_design) %>%
#     mutate(across(starts_with("y"), \(y) as.numeric(y) - 1))
#   ystart <- min(grep("^y[0-9]", names(v$variables))) - 1
#   idx <- combn(p, 2)
#   for (k in seq_len(ncol(idx))) {
#     i <- idx[1, k]
#     j <- idx[2, k]
#     varname <- paste0("y", i, ".", j, collapse = "")
#     yi <- v$variables[, i + ystart, drop = TRUE]
#     yj <- v$variables[, j + ystart, drop = TRUE]
#     yij <- (yi == 1) * (yj == 1)  # both positive
#     v$variables[[varname]] <- yij
#     y_form <- paste0(y_form, paste0(" + ", varname))
#   }
#
#   x <- v$variables[, -(1:ystart)]
#   survey::svyvar(x, .svy_design) %>%
#     as.matrix()  # NOT GOOD BECAUSE svyvar() not as efficient as cov.wt()
#   # wt <- 1 / .svy_design$prob
#   # res <- cov.wt(x, wt, center = c(pidot1, pidot2))$cov
#   # res
# }

# create_Sigma_univariate_matrix_complex <- function(.lavobject, .svy_design) {
#   list2env(extract_lavaan_info(.lavobject), environment())
#   list2env(get_uni_bi_moments(.lavobject), environment())
#
#   # Get strata information
#   strata  <- .svy_design$strata[, 1, drop = FALSE]
#   colnames(strata) <- "strata"
#   a <- unique(strata[, 1])
#
#   # Get cluster information
#   cluster <- .svy_design$cluster
#   colnames(cluster)[1] <- "cluster"
#
#   dat <-
#     dat %>%
#     bind_cols(strata, cluster)
#
#   # Create the data set of univariate and bivariate responses (1/0)
#   tmp <- dat %>%
#     select(starts_with("y")) %>%
#     mutate(across(starts_with("y"), function(x) as.numeric(x == 2)))
#   tmp <- bind_cols(dat %>% select(-starts_with("y")), tmp)
#
#   Sigma <- list()
#   for (i in a) {
#     # Filter the strata a
#     current_dat <- tmp %>%
#       filter(strata == i)
#
#     # Obtain the responses in this strata (i.e. according to S_ab)
#     y <- current_dat %>%
#       select(starts_with("y")) %>%
#       mutate(across(starts_with("y"), function(x) as.numeric(x == 1))) %>%
#       as.matrix()
#
#     uab <- {sweep(y, 2, pidot1, "-") * current_dat$wt} %>%
#       as_tibble() %>%
#       bind_cols(current_dat %>% select(-starts_with("y")), .) %>%
#       group_by(cluster) %>%
#       summarise(across(starts_with("y"), sum), .groups = "drop") %>%
#       select(starts_with("y")) %>%
#       as.matrix()
#     na <- nrow(uab)
#
#     ubar <- apply(uab, 2, mean)
#
#     Sigma[[i]] <- crossprod(sweep(uab, 2, ubar, "-")) * na / (na - 1)
#   }
#   Reduce("+", Sigma) / N
# }

get_w_for_delta <- function(.lavobject, .svy_design) {
  if (is.null(.svy_design)) return(NA)

  list2env(extract_lavaan_info(.lavobject), environment())
  list2env(get_uni_bi_moments(.lavobject), environment())
  y_form <- paste0("~ ", paste0("y", 1:p, collapse = " + "), sep = "")
  v <-
    srvyr::as_survey(.svy_design) %>%
    mutate(across(starts_with("y"), \(y) as.numeric(y) - 1))
  ystart <- min(grep("^y[0-9]", names(v$variables))) - 1
  idx <- combn(p, 2)
  for (k in seq_len(ncol(idx))) {
    i <- idx[1, k]
    j <- idx[2, k]
    varname <- paste0("y", i, ".", j, collapse = "")
    yi <- v$variables[, i + ystart, drop = TRUE]
    yj <- v$variables[, j + ystart, drop = TRUE]
    yij <- (yi == 1) * (yj == 1)  # both positive
    v$variables[[varname]] <- yij
    y_form <- paste0(y_form, paste0(" + ", varname))
  }
  x <- v$variables[, -(1:ystart)]
  wt <- 1 / .svy_design$prob
  apply(x, 2, \(x) sum(x * wt ^ 2)) / apply(x, 2, \(x) sum(x * wt))
}

## ---- Test preliminaries ----------------------------------------------------
calc_test_stuff <- function(lavobject, svy_design = NULL, .H_inv = NULL,
                            .Delta_mat_list  = NULL, .Sigma2  = NULL, .wstar = 1,
                            .pi_tilde = NULL, bootstrap = FALSE, nboot = 100) {
  # These are all the stuff needed to compute the test statistic X2
  list2env(extract_lavaan_info(lavobject), environment())

  if (is.null(.H_inv)) {
    H_inv <- get_sensitivity_inv_mat(.lavobject = lavobject)
  } else {
    H_inv <- .H_inv
  }
  if (is.null(.Delta_mat_list)) {
    Delta_mat_list <- get_Delta_mats(.lavobject = lavobject)
  } else {
    Delta_mat_list <- .Delta_mat_list
  }
  if (is.null(.Sigma2)) {
    if (is.null(svy_design)) {
      Sigma2 <- create_Sigma2_matrix(lavobject)
    } else {
      Sigma2 <- create_Sigma2_matrix_complex(lavobject, svy_design, bootstrap,
                                             nboot)
    }
  } else {
    Sigma2 <- .Sigma2
  }
  if (is.null(.pi_tilde)) {
    pi_tilde <- unlist(lav_tables_pairwise_model_pi(lavobject))
  } else {
    pi_tilde <- .pi_tilde
  }

  wstar_tmp <- attr(svy_design$variables, "wstar")
  if (!is.null(wstar_tmp)) {
    wstar <- wstar_tmp
    wstar <- 1  # Actually I think this should just be 1...
  } else {
    wstar <- .wstar
  }

  # Delta tilde and Delta_2 matrices
  Delta2 <- Delta_mat_list$Delta2
  Delta_til <- Delta_mat_list$Delta_til

  # Joint proportions, model probabilities and Sigma matrices
  list2env(get_uni_bi_moments(lavobject), environment())
  p2_hat <- c(pdot1, pdot2)
  pi2_hat <- c(pidot1, pidot2)

  # Residuals
  e2_hat <- p2_hat - pi2_hat

  # B matrix converting parameters to residuals
  B_mat <- Beta_mat_design(p)  # Myrsini's B_mat is actually G %*% ginv(T2)
  # Bivariate pairs are usually ordered as follows: y1y2, y1y3, ..., y1yp,
  # y2y3, y2y4, ..., y2yp, ..., yp-1yp. In each pair, the permutation is 00,
  # 10, 01, 11. Thus, every four entry sums to 1.
  B2 <- wstar * H_inv %*% t(Delta_til * (1 / pi_tilde)) %*% B_mat

  # Omega_2 matrix
  Omega2 <- Sigma2 -
    Delta2 %*% B2 %*% Sigma2 -
    Sigma2 %*% t(B2) %*% t(Delta2) +
    Delta2 %*% B2 %*% Sigma2 %*% t(B2) %*% t(Delta2)
  Omega2 <- (Omega2 + t(Omega2)) / 2

  # "Approximate" Omega_2 matrix using Fisher information
  Fisher_mat_inv <- N * get_sensitivity_inv_mat(lavobject, "Sandwich")
  Omega2_approx <- Sigma2 - Delta2 %*% Fisher_mat_inv %*% t(Delta2)
  Omega2_approx <- (Omega2_approx + t(Omega2_approx)) / 2

  res <- list(
    N = N,               # Sample size
    q = ncol(H_inv),     # No. of parameters
    p = p,               # No. of items
    S = length(e2_hat),  # p(p+1)/2

    # Residuals and the covariance matrix of the residuals
    e2_hat        = e2_hat,
    Omega2        = Omega2,
    Omega2_approx = Omega2_approx,

    # Proportions and model probabilities
    p2_hat  = p2_hat,
    pi2_hat = pi2_hat,

    # Other matrices
    Sigma2    = Sigma2,
    H_inv     = H_inv,
    B2        = B2,
    Delta2    = Delta2,

    # Asparouhov Muthen approx delta
    approx_delta = get_w_for_delta(lavobject, svy_design)
  )
  attr(res, "bingof_calc_test_stuff") <- TRUE
  res
}

test_begin <- function(.lavobject, .approx_Omega2, .svy_design = NULL) {
  if (is(.lavobject, "lavaan")) {
    the_test_stuff <- calc_test_stuff(.lavobject, .svy_design)
  } else {
    the_test_stuff <- .lavobject
  }

  # Should we use approximate Omega2, which is the Godambe sandwich estimator
  # (inverse) for the Fisher information instead of sensitivity matrix?
  if (isTRUE(.approx_Omega2)) {
    Covmat <- the_test_stuff$Omega2_approx
  } else {
    Covmat <- the_test_stuff$Omega2
  }
  the_test_stuff$Omega2 <- Covmat

  the_test_stuff
}

after_test <- function(x, .Xi, .S) {
  # x is a data.frame containing X2 and df
  x %>%
    mutate(pval = pchisq(X2, df, lower.tail = FALSE),
           Xi_rank = Matrix::rankMatrix(.Xi),
           S = .S)
}

moment_match <- function(X2, Xi, Omega2, df = NULL, order) {
  XiOmega2 <- Xi %*% Omega2
  mu1 <- sum(diag(XiOmega2))
  mu2 <- 2 * sum(diag(XiOmega2 %*% XiOmega2))
  mu3 <- 8 * sum(diag(XiOmega2 %*% XiOmega2 %*% XiOmega2))

  order <- match.arg(as.character(order), c("0", "1", "2", "3"))
  if (order == "3") {
    b <- mu3  / (4 * mu2)
    c <- mu2 / (2 * b ^ 2)
    a <- mu1 - b * c
  }
  # if (a > X2) order <- "2"
  if (order == "2") {
    b <- mu2 / (2 * mu1)
    c <- mu1 / b
    a <- 0
  }
  if (order == "1" & !is.null(df)) {
    c <- df
    b <- mu1 / c
    a <- 0
  }
  if (order == "0" & !is.null(df)) {
    c <- df
    b <- 1
    a <- 0
    cli::cli_warn("Test statistics should be moment matched. This option is provided for testing purposes only.")
  }

  data.frame(X2 = (X2 - a) / b, df = c)
}

#' Limited information goodness-of-fit tests
#'
#' @name ligof-test-stats
#' @rdname ligof-test-stats
#'
#' @param object A [lavaan::lavaan()] fit object.
#' @param approx_Omega2 `r lifecycle::badge('experimental')` (logical)  Should an approximate
#'   residual covariance matrix \eqn{\Omega_2} be used? Defaults to `FALSE`.
#' @param svy_design (optional) A [survey::svydesign()] object.
#' @param .order (integer) Either the number of moments to match for the
#'   chi-square test statistic matching procedure (choose from 1--3), or the
#'   Rao-Scott type adjustment order (choose from 1 or 2).
#'
#' @returns A data frame containing the test statistics \eqn{X^2}, degrees of
#'   freedom, name of the test, and the \eqn{p}-value.
#' @seealso [all_tests()]
#'
#' @examples
#' fit <- lavaan::sem(txt_mod(1), gen_data_bin(1, n = 500), std.lv = TRUE,
#'                    estimator = "PML")
#' Wald_test(fit)
NULL

#' @describeIn ligof-test-stats The Wald test statistic.
#' @export
Wald_test <- function(object, approx_Omega2 = FALSE, svy_design = NULL) {
  list2env(test_begin(object, approx_Omega2, svy_design), environment())

  Xi <- MASS::ginv(Omega2)
  X2 <- N * colSums(e2_hat * (Xi %*% e2_hat))
  data.frame(X2 = X2, df = Matrix::rankMatrix(Omega2) - q, name = "Wald") %>%
    after_test(., Xi, S)
}

Wald_test_v2 <- function(object, approx_Omega2 = FALSE, svy_design = NULL,
                         .order = 3) {
  list2env(test_begin(object, approx_Omega2, svy_design), environment())

  Xi <- diag(1 / diag(Omega2))
  # X2 <- N * colSums(e2_hat * (Xi %*% e2_hat))
  X2 <- N * sum(e2_hat / diag(Omega2) * e2_hat)
  out <- moment_match(X2, Xi, Omega2, df = S - q, order = .order)
  cbind(out, name = paste0("WaldDiag,MM", .order)) %>%
    after_test(., Xi, S)
}

#' @describeIn ligof-test-stats The Wald test statistic using a simple diagonal
#'   \eqn{\Omega_2} matrix.
#' @export
Wald_diag_test <- Wald_test_v2

Wald_diag_RS_test <- function(object, approx_Omega2 = FALSE, svy_design = NULL,
                              .order = 2) {
  .order <- match.arg(as.character(.order), c("1", "2"))
  list2env(test_begin(object, approx_Omega2, svy_design), environment())

  omega_diag <- diag(Omega2)
  Xi <- diag(1 / omega_diag)
  X2 <- N * colSums(e2_hat * (Xi %*% e2_hat))

  # Rao-Scott adjustment
  mat <- Omega2 %*% diag(1 / omega_diag)
  delta <- eigen(mat)$values

  X2 <- X2 / mean(delta)
  df <- S
  if (.order == "2") {
    a_sq <- mean((delta - mean(delta)) ^ 2) / mean(delta) ^ 2
    # a_sq <- (sd(delta) / mean(delta)) ^ 2
    X2 <- X2 / (1 + a_sq)
    df <- S / (1 + a_sq)
  }

  data.frame(X2 = X2, df = df, name = paste0("WaldDiag,RS", .order)) %>%
    after_test(., Xi, S)
}

Wald_diag_RS_approx_test <- function(object, approx_Omega2 = FALSE,
                                     svy_design = NULL, .order = 2) {
  .order <- match.arg(as.character(.order), c("1", "2"))
  list2env(test_begin(object, approx_Omega2, svy_design), environment())

  omega_diag <- diag(Omega2)
  Xi <- diag(1 / omega_diag)
  X2 <- N * colSums(e2_hat * (Xi %*% e2_hat))
  delta <- approx_delta

  X2 <- X2 / mean(delta)
  df <- S
  if (.order == "2") {
    a_sq <- mean((delta - mean(delta)) ^ 2) / mean(delta) ^ 2
    # a_sq <- (sd(delta) / mean(delta)) ^ 2
    X2 <- X2 / (1 + a_sq)
    df <- S / (1 + a_sq)
  }

  data.frame(X2 = X2, df = df, name = paste0("WaldDiag,RSaprx", .order)) %>%
    after_test(., Xi, S)
}

Wald_test_v3 <- function(object, svy_design = NULL) {
  list2env(test_begin(object, .approx_Omega2 = FALSE, svy_design),
           environment())

  Delta2comp <- mcompanion::null_complement(Delta2)
  Xi <- Delta2comp %*% MASS::ginv(
    t(Delta2comp) %*% Sigma2 %*% Delta2comp
  ) %*% t(Delta2comp)
  X2 <- N * colSums(e2_hat * (Xi %*% e2_hat))
  data.frame(X2 = X2, df = S - q, name = "WaldVCF") %>%
    after_test(., Xi, S)
}

#' @describeIn ligof-test-stats The Wald test statistic bypassing the
#'   \eqn{\Omega_2} matrix (uses orthogonal complements of \eqn{\Delta_2}).
#' @export
Wald_vcovf_test <- Wald_test_v3

Pearson_test_v1 <- function(object, approx_Omega2 = FALSE, svy_design = NULL,
                            .order = 2) {
  .order <- match.arg(as.character(.order), c("1", "2"))
  list2env(test_begin(object, approx_Omega2, svy_design), environment())

  Xi <- diag(1 / pi2_hat)
  X2 <- N * colSums(e2_hat * (Xi %*% e2_hat))

  # Rao-Scott adjustment
  mat <- Omega2 %*% diag(1 / pi2_hat)
  delta <- eigen(mat)$values

  X2 <- X2 / mean(delta)
  df <- S
  if (.order == "2") {
    a_sq <- mean((delta - mean(delta)) ^ 2) / mean(delta) ^ 2
    # a_sq <- (sd(delta) / mean(delta)) ^ 2
    X2 <- X2 / (1 + a_sq)
    df <- S / (1 + a_sq)
  }

  data.frame(X2 = X2, df = df, name = paste0("Pearson,RS", .order)) %>%
    after_test(., Xi, S)
}

#' @describeIn ligof-test-stats The Pearson test with \eqn{p}-values calculated using a Rao-Scott type adjustment.
#' @export
Pearson_RS_test <- Pearson_test_v1

Pearson_RS_approx_test <- function(object, approx_Omega2 = FALSE,
                                   svy_design = NULL, .order = 2) {
  .order <- match.arg(as.character(.order), c("1", "2"))
  list2env(test_begin(object, approx_Omega2, svy_design), environment())

  Xi <- diag(1 / pi2_hat)
  X2 <- N * colSums(e2_hat * (Xi %*% e2_hat))
  delta <- approx_delta

  X2 <- X2 / mean(delta)
  df <- S
  if (.order == "2") {
    a_sq <- mean((delta - mean(delta)) ^ 2) / mean(delta) ^ 2
    # a_sq <- (sd(delta) / mean(delta)) ^ 2
    X2 <- X2 / (1 + a_sq)
    df <- S / (1 + a_sq)
  }

  data.frame(X2 = X2, df = df, name = paste0("Pearson,RSaprx", .order)) %>%
    after_test(., Xi, S)
}

# Pearson_test_v3 <- function(object, approx_Omega2 = FALSE, svy_design = NULL,
#                             .order = 2) {
#   .order <- match.arg(as.character(.order), c("1", "2"))
#   list2env(test_begin(object, approx_Omega2, svy_design), environment())
#
#   Xi <- diag(1 / pi2_hat)
#   X2 <- N * colSums(e2_hat * (Xi %*% e2_hat))
#
#   # Rao-Scott adjustment
#   tmp <- eigen(Omega2)
#   U <- diag(tmp$val)
#   V <- tmp$vec
#   Omegahalf <- V %*% sqrt(U) %*% t(V)
#   # same as below
#   Omegahalf <- t(chol(Omega2))
#   mat <- t(Omegahalf) %*% Xi %*% (Omegahalf)
#   delta <- eigen(mat)$values
#
#   X2 <- X2 / mean(delta)
#   df <- S
#   if (.order == "2") {
#     a_sq <- mean((delta - mean(delta)) ^ 2) / mean(delta) ^ 2
#     # a_sq <- (sd(delta) / mean(delta)) ^ 2
#     X2 <- X2 / (1 + a_sq)
#     df <- S / (1 + a_sq)
#   }
#
#   data.frame(X2 = X2, df = df, name = "PearsonRS") %>%
#     after_test(., Xi, S)
# }

Pearson_test_v2 <- function(object, approx_Omega2 = FALSE, svy_design = NULL,
                            .order = "3") {
  list2env(test_begin(object, approx_Omega2, svy_design), environment())

  Xi <- diag(1 / pi2_hat)
  X2 <- N * colSums(e2_hat * (Xi %*% e2_hat))
  out <- moment_match(X2, Xi, Omega2, df = S - q, order = .order)
  cbind(out, name = paste0("Pearson,MM", .order)) %>%
    after_test(., Xi, S)
}

#' @describeIn ligof-test-stats The Pearson test with \eqn{p}-values calculated using a moment-matching procedure.
#' @export
Pearson_test <- Pearson_test_v2

#' @describeIn ligof-test-stats The residual sum of squares (RSS) test. Uses moment-matching for \eqn{p}-value calculations.
#' @export
RSS_test <- function(object, approx_Omega2 = FALSE, svy_design = NULL,
                     .order = "3") {
  list2env(test_begin(object, approx_Omega2, svy_design), environment())

  Xi <- diag(S)
  # X2 <- N * colSums(e2_hat * (Xi %*% e2_hat))
  X2 <- N * sum(e2_hat ^ 2)
  out <- moment_match(X2, Xi, Omega2, df = S - q, order = .order)
  cbind(out, name = paste0("RSS,MM", .order)) %>%
    after_test(., Xi, S)
}

#' @describeIn ligof-test-stats The multinomial test. Uses moment-matching for \eqn{p}-value calculations.
#' @export
Multn_test <- function(object, approx_Omega2 = FALSE, svy_design = NULL,
                       .order = "3") {
  list2env(test_begin(object, approx_Omega2, svy_design), environment())

  Xi <- MASS::ginv(Sigma2)
  X2 <- N * colSums(e2_hat * (Xi %*% e2_hat))
  out <- moment_match(X2, Xi, Omega2, df = S - q, order = .order)
  cbind(out, name = paste0("Multn,MM", .order)) %>%
    after_test(., Xi, S)
}

#' Return all test statistics values
#'
#' @inherit ligof-test-stats params return
#' @param sim (integer) Optional and used for large-scale simulations.
#' @param bootstrap (boolean) Use the `svyrep` package to compute bootstrap
#'   variance estimator?
#' @param nboot (integer) Optional if `bootstrap = TRUE` then how many bootstrap
#'   replications?
#' @param Sigma2 (for internal testing only)
#' @param Hinv (for internal testing only)
#' @param wstar (for internal testing only)
#'
#' @returns Additionally, if `sim` argument is provided, two columns are
#'   appended: Whether the [lavaan::lavaan()]  fit has `converged` and the
#'   matrix rank of \eqn{\Omega_2} (useful to see if any computational issues
#'   arose during model fit.)
#' @export
#'
#' @examples
#' fit <- lavaan::sem(txt_mod(1), gen_data_bin(1, n = 500), std.lv = TRUE,
#'                    estimator = "PML")
#' all_tests(fit)
all_tests <- function(object, svy_design = NULL, sim = NULL, Sigma2 = NULL,
                      wstar = 1, bootstrap = FALSE, nboot = 100, Hinv = NULL) {
  if (isTRUE(attr(object, "bingof_calc_test_stuff"))) {
    test_stuff <- object
  } else {
    test_stuff <- calc_test_stuff(object, svy_design, .Sigma2 = Sigma2,
                                  .wstar = wstar, bootstrap = bootstrap,
                                  nboot = nboot, .H_inv = Hinv)
  }

  res <- bind_rows(
    Wald_test(test_stuff),
    Wald_vcovf_test(test_stuff),
    # Wald_diag_test(test_stuff, .order = 1),
    # Wald_diag_test(test_stuff, .order = 2),
    Wald_diag_test(test_stuff, .order = 3),
    # Wald_diag_RS_test(test_stuff, .order = 1),
    Wald_diag_RS_test(test_stuff, .order = 2),
    # Wald_diag_RS_approx_test(test_stuff, .order = 2),

    # Pearson_test(test_stuff, .order = 1),
    # Pearson_test(test_stuff, .order = 2),
    Pearson_test(test_stuff, .order = 3),
    # Pearson_RS_test(test_stuff, .order = 1),
    Pearson_RS_test(test_stuff, .order = 2)
    # Pearson_RS_approx_test(test_stuff, .order = 2)

    # RSS_test(test_stuff, .order = 1),
    # RSS_test(test_stuff, .order = 2),
    # RSS_test(test_stuff, .order = 3),
    # Multn_test(test_stuff, .order = 1),
    # Multn_test(test_stuff, .order = 2),
    # Multn_test(test_stuff, .order = 3)
  ) %>%
    as_tibble()
  if (!is.null(sim)) {
    res <- bind_cols(tibble(i = sim), res) %>%
      bind_cols(
        converged = object@Fit@converged,
        Omega2_rank = Matrix::rankMatrix(test_stuff$Omega2)
      ) %>%
      suppressWarnings()
  }

  res
}
