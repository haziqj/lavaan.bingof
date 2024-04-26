# res <- fit_facmod_pml(1, samp = "clust")
create_Sigma2_matrix_complex <- function(.lavobject, strat, clust) {
  list2env(extract_lavaan_info(.lavobject), environment())
  list2env(get_uni_bi_moments(.lavobject), environment())
  p2_hat <- c(pdot1, pdot2)     # the sample proportions?
  pi2_hat <- c(pidot1, pidot2)  # or uni and bivariate moments (model implied)?

  method <- "strcl"
  if (missing(strat)) method <- "clust"
  if (missing(clust)) method <- "strat"

  # Prepare the sample data --------------------------------------------------
  dat <-
    dat %>%
    mutate(across(everything(), \(y) as.numeric(y) - 1))
  idx <- combn(p, 2)
  varnames <- colnames(dat)
  for (k in seq_len(ncol(idx))) {
    # Create the uni and bivariate data
    i <- idx[1, k]
    j <- idx[2, k]
    varname <- paste0(varnames[i], varnames[j])
    yi <- dat[, i, drop = TRUE]
    yj <- dat[, j, drop = TRUE]
    yij <- (yi == 1) * (yj == 1)  # both positive
    dat[[varname]] <- yij
  }
  dat <- scale(dat, center = pi2_hat, scale = FALSE)

  if (method == "strat") {
    strata <- unique(strat)
    nstrat <- length(strata)
    za <- list()
    for (a in seq_along(strata)) {
      dat_a <- dat[strat == strata[a], ]
      wt_a <- wt[strat == strata[a]]
      tmp <- apply(dat_a * wt_a, 2, sum) / sum(wt_a)
      za[[a]] <- tcrossprod(tmp) * nrow(dat_a) / (nrow(dat_a) - 1)
    }
    Sigma2 <- Reduce(`+`, za)
    print("METHOD = STRAT")
  }

  if (method == "clust") {
    clusters <- unique(clust)
    nclust <- length(clusters)
    zb <- list()
    for (b in seq_along(clusters)) {
      dat_b <- dat[clust == clusters[b], ]
      wt_b <- wt[clust == clusters[b]]

      zb[[b]]  <- apply(dat_b * wt_b, 2, sum) / sum(wt_b)
    }
    zbar <- apply(do.call(cbind, zb), 1, mean)
    ZZt <-
      lapply(zb, \(z) tcrossprod(z - zbar)) |>
      Reduce(f = `+`)
    Sigma2 <- nclust / (nclust - 1) * ZZt
    print("METHOD = CLUST")
  }

  if (method == "strcl") {
    strata <- unique(strat)
    nstrat <- length(strata)
    za <- list()
    for (a in seq_along(strata)) {
      dat_a <- dat[strat == strata[a], ]
      wt_a <- wt[strat == strata[a]]
      clust_a <- clust[strat == strata[a]]
      clusters_a <- unique(clust_a)
      nclust_a <- length(clusters_a)
      zab <- list()
      for (b in seq_along(clusters_a)) {
        dat_ab <- dat_a[clust_a == clusters_a[b], ]
        wt_ab <- wt_a[clust_a == clusters_a[b]]
        zab[[b]] <- apply(dat_ab * wt_ab, 2, sum) / sum(wt_ab)
      }
      zbar_a <- apply(do.call(cbind, zab), 1, mean)
      ZZt <-
        lapply(zab, \(z) tcrossprod(z - zbar_a)) |>
        Reduce(f = `+`)
      za[[a]] <- nclust_a / (nclust_a - 1) * ZZt
    }
    Sigma2 <- Reduce(`+`, za)
    print("METHOD = STRCL")
  }

  return(Sigma2)
}


# Sigma2_pop   <- attr(make_population(1, Sigma2_attr = TRUE), "Sigma2")
# Sigma2_theor <- create_Sigma2_matrix(res$fit, method = "theoretical")
# Sigma2_clust <- create_Sigma2_matrix_complex(res$fit, strat = res$dat$type, clust = res$dat$school)
