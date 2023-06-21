# population should have 50 strata of unequal sizes Nh. Total pop size is roughly N = 1,000,000.
# strat sample from each strata, select say 20 individuals. The prob of selection is 20 / Nh

make_population2 <- function(model_no = 1, seed = 123, H1 = FALSE,
                             return_all = FALSE) {

  set.seed(seed)

  # Population settings --------------------------------------------------------
  nschools       <- c(400, 1000, 600)  # types A, B, C
  avg_class_size <- c(15, 25, 20)
  avg_nstudent   <- 500  # so pop total is nschools x avg_nstudent
  sd_nstudent    <- 100

  # Set up the loadings and covariance matrices --------------------------------
  Lambda      <- loading_mat(model_no)
  neta        <- ncol(Lambda)  # q
  nitems      <- nrow(Lambda)  # p
  Psi         <- cov_lv_mat(model_no)
  Theta       <- matrix(0, nrow = nitems, ncol = nitems)
  diag(Theta) <- 1 - diag(Lambda %*% Psi %*% t(Lambda))
  tau         <- get_tau(model_no)

  # Prepare population data ----------------------------------------------------
  letters <- c(letters, sapply(letters, function(x) paste0(x, letters)))
  pop <-
    tibble(type = LETTERS[1:3],
           nschools = nschools,
           avg_class_size = avg_class_size) %>%
    rowwise() %>%
    # mutate(school = list(sprintf("%04d", seq_len(nschools)))) %>%
    mutate(school = list(seq_len(nschools))) %>%
    unnest_longer(school) %>%
    mutate(nstudents = round(rnorm(dplyr::n(), avg_nstudent,
                                   sd = sd_nstudent))) %>%
    rowwise() %>%
    mutate(class = list(sample(letters[seq_len(nstudents / avg_class_size)],
                               size = nstudents, replace = TRUE))) %>%
    unnest_longer(class) %>%
    select(type, school, class) %>%
    arrange(type, school, class) %>%
    mutate(school = paste0(type, school),
           class = paste0(school, class))
  N <- nrow(pop)

  # Generate the data ----------------------------------------------------------
  eta     <- mvnfast::rmvn(n = N, mu = rep(0, neta), sigma = Psi)
  epsilon <- mvnfast::rmvn(n = N, mu = rep(0, nitems), sigma = Theta)
  # eta     <- mnormt::rmnorm(n = N, mean = rep(0, neta), varcov = Psi)
  # epsilon <- mnormt::rmnorm(n = N, mean = rep(0, nitems), varcov = Theta)
  ystar   <- tcrossprod(eta, Lambda) + epsilon

  if (isTRUE(H1)) {
    # Add an extra factor to misspecify the model fit (for power simulations)
    extra_Lambda <- Lambda[, 1, drop = FALSE] + rnorm(nitems, sd = 0.1)
    if (model_no <= 3) { extra_Lambda[seq(2, nitems, by = 2), 1] <- 0 }
    ystar <- ystar + t(extra_Lambda %*% rnorm(N))
    ystar <- scale(ystar)
  }

  # repair eta
  eta <- as.data.frame(eta)
  colnames(eta) <- paste0("eta", seq_len(neta))

  # Stratify according to latent variable --------------------------------------
  abil_order <-
    as_tibble(eta) %>%
    mutate(z = rowSums(across(everything())),
           rn = row_number()) %>%
    select(z, rn) %>%
    arrange(desc(z)) %>%
    mutate(type = pop$type) %>%
    group_by(type) %>%
    mutate(rn2 = sample(dplyr::n(), replace = FALSE)) %>%
    arrange(type, rn2) %>%
    pull(rn)
  ystar <- ystar[abil_order, ]

  # Get the response patterns --------------------------------------------------
  y <-
    t(apply(ystar, 1, function(x) as.numeric(x > tau))) %>%
    as.data.frame()
  colnames(y) <- paste0("y", seq_len(nitems))
  Sigma2 <- (N - 1) / N * cov(convert_dat_to_unibiv(y))

  ystar <- as.data.frame(ystar)
  colnames(ystar) <- paste0("ystar", seq_len(nitems))

  if (isTRUE(return_all)) {
    res <- bind_cols(pop, y, ystar, eta[abil_order, , drop = FALSE])
  } else {
    res <- bind_cols(pop, y)
  }

  attr(res, "Sigma2") <- Sigma2
  res
}


#' @export
gen_data_bin_strat2 <- function(population = make_population(1, seed = NULL),
                                  npsu = 2, n, seed = NULL) {
  # I want to sample x number of schools, and select all students from that
  # school.

  # 1-stage stratified sampling
  set.seed(seed)
  if (!missing(n)) npsu <- max(1, round(n / 1500, 0))

  # Weights
  school_info <-
    population %>%
    group_by(type) %>%
    summarise(nschools = n_distinct(school)) %>%
    mutate(
      prob = npsu / nschools,
      wt = 1 /prob
    )

  # Sampling of PSUs
  psu_sampled <- population %>%
    group_by(type) %>%
    distinct(school) %>%
    slice_sample(n = 1) %>%
    ungroup()

  # Add the classrooms to the sample
  sampled <-
    inner_join(population, psu_sampled, by = c("type", "school")) %>%
    left_join(school_info, by = c("type")) %>%
    select(type, school, class, wt, starts_with("y")) %>%
    mutate(wt = wt / sum(wt) * dplyr::n(),
           across(starts_with("y"), ordered)) %>%
    arrange(type, school, class)

  sampled
}
