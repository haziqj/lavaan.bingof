convert_dat_to_unibiv <- function(dat) {
  p <- ncol(dat)
  idx <- combn(p, 2)

  for (k in seq_len(ncol(idx))) {
    i <- idx[1, k]
    j <- idx[2, k]
    varname <- paste0("y", i, ".", j, collapse = "")
    yi <- dat[, i, drop = TRUE]
    yj <- dat[, j, drop = TRUE]
    yij <- (yi == 1) * (yj == 1)  # both positive
    dat[[varname]] <- yij
  }
  dat
}

#' Simulate the school population data set
#'
#' @inheritParams gen_data_bin
#' @param return_all (logical) Return the underlying latent variables (\eqn{y^*}
#'   and \eqn{\eta}) as well? well?
#' @param Sigma2_attr (logical) Should the population Sigma2 matrix be
#'   computed and stored as an attribute?
#'
#' @return A [tibble()] containing ordinal binary values (0/1) for the items, as
#'   well as the population stratum and clusters (`type`, `school`, `class`).
#' @export
#'
#' @seealso [gen_data_bin_strat()], [gen_data_bin_clust()],
#'   [gen_data_bin_strcl()], [gen_data_bin_srs()]
#'
#'
#' @examples
#' \dontrun{
#' make_population(1)
#' }
#'
make_population <- function(model_no = 1, seed = 123, H1 = FALSE,
                            return_all = FALSE, Sigma2_attr = FALSE) {

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
    if (model_no == 1) { extra_Lambda[seq(2, nitems, by = 2), 1] <- 0 }
    if (model_no == 2) { extra_Lambda[seq(2, nitems, by = 4), 1] <- 0 }
    if (model_no == 3) { extra_Lambda[seq(2, nitems, by = 6), 1] <- 0 }
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
    #   # t(apply(ystar, 1, function(x) as.numeric(x > tau))) %>%
    1 * (ystar > matrix(tau, nrow = N, ncol = nitems, byrow = TRUE)) %>%
    as.data.frame()
  colnames(y) <- paste0("y", seq_len(nitems))

  ystar <- as.data.frame(ystar)
  colnames(ystar) <- paste0("ystar", seq_len(nitems))

  if (isTRUE(return_all)) {
    res <- bind_cols(pop, y, ystar, eta[abil_order, , drop = FALSE])
  } else {
    res <- bind_cols(pop, y)
  }

  Sigma2 <- NULL
  if (isTRUE(Sigma2_attr)) {
    Sigma2 <- (N - 1) / N * cov(convert_dat_to_unibiv(y))
  }
  attr(res, "Sigma2") <- Sigma2

  res
}

#' Sample from the school population
#'
#' @description There are several functions to simulate a complex (or even
#'   simple random) sampling procedure on the school population generated by
#'   [make_population()]:
#'
#' - `gen_data_bin_srs()` performs a **simple random sample** on the entire population.
#' - `gen_data_bin_strat()` performs a **stratified** sampling procedure with the school `type` as strata, and the students as the PSU.
#'
#' - `gen_data_bin_clust()` performs a **two-stage cluster** sampling procedure (ignoring the strata) with `school`s being the PSU, and further clustering of `class`.
#'
#' - `gen_data_bin_strcl()` performs a **two-stage stratified cluster** sampling procedure where students are nested within `class` within `school` (PSU) within `type` (strata).
#'
#' @inherit make_population params return
#' @param population (tibble) Population generated by [make_population()].
#' @param npsu (integer > 0) For the stratified sampling, this is the number of
#'   students within each strata to sample. For the cluster sampling, this is
#'   the number of schools (clusters) to sample. For the stratified cluster
#'   sampling, this is also the number of schools (clusters) per school type
#'   (strata). For the SRS procedure, this is the exact sample size.
#' @param n (optional,integer > 0) Sample size. If provided, then the `npsu`
#'   argument is ignored and adjusted accordingly to achieve a sample size of
#'   (approximately) `n`.
#'
#' @seealso [make_population()]
#'
#' @name gen_data_bin_complex
#' @rdname gen_data_bin_complex
#'
#' @examples
#' \dontrun{
#' pop <- make_population(2)
#' gen_data_bin_srs(pop)
#' gen_data_bin_strat(pop)
#' gen_data_bin_clust(pop)
#' gen_data_bin_strcl(pop)
#' }
#'
#'
NULL

#' @rdname gen_data_bin_complex
#' @export
gen_data_bin_srs <- function(population = make_population(1, seed = NULL),
                             npsu = 3000, n, seed = NULL) {
  set.seed(seed)
  if (!missing(n)) npsu <- n
  slice_sample(population, n = npsu, replace = FALSE) %>%
    mutate(across(starts_with("y"), ordered))
}

gen_data_bin_complex1 <- function(population = make_population(1, seed = NULL),
                                  npsu = 1000, n, seed = NULL) {
  # 1-stage stratified sampling
  set.seed(seed)
  if (!missing(n)) npsu <- round(n / 3, 0)

  # Weights
  school_info <-
    population %>%
    group_by(type) %>%
    summarise(
      students_in_school_type = dplyr::n(),
    ) %>%
    mutate(
      prob = npsu / students_in_school_type,
      wt = 1 / prob
    )

  # Add back the school_info to the population
  population <-
    population %>%
    left_join(school_info, by = c("type"))

  # Sampling of PSUs
  sampled <- population %>%
    group_by(type) %>%
    slice_sample(n = npsu, replace = FALSE) %>%
    ungroup() %>%
    select(-starts_with("ystar")) %>%
    select(type, school, class, wt, starts_with("y")) %>%
    mutate(wt = wt / sum(wt) * dplyr::n(),
           across(starts_with("y"), ordered)) %>%
    arrange(type, school, class)

  sampled
}

#' @rdname gen_data_bin_complex
#' @export
gen_data_bin_strat <- gen_data_bin_complex1

gen_data_bin_complex2 <- function(population = make_population(1, seed = NULL),
                                  npsu = 140, n, seed = NULL) {
  # 2-stage cluster sampling
  set.seed(seed)
  if (!missing(n)) npsu <- round(n / 21.5, 0)

  # Sampling of PSUs
  pop2 <-
    population %>%
    select(type, school, class) %>%
    group_by(type, school) %>%
    summarise(
      nstudents = dplyr::n(),
      pr_class_selected = 1 / length(unique(class)),
      .groups = "drop"
    ) %>%
    mutate(
      pr_school_selected = nstudents / sum(nstudents),
      prob = pr_school_selected * pr_class_selected,
      wt = 1 / prob
    )
  # mutate(pr_school_selected = npsu / dplyr::n()) %>%

  psu_sampled <-
    pop2 %>%
    slice_sample(n = npsu, weight_by = nstudents, replace = FALSE) %>%
    arrange(type, school)

  # Sampling of classes within PSUs (using SRS)
  classes_sampled <-
    inner_join(population, psu_sampled, by = c("type", "school")) %>%
    distinct(type, school, class) %>%
    group_by(type, school) %>%
    slice_sample(n = 1, replace = FALSE) %>%
    ungroup()

  sampled <-
    inner_join(population, classes_sampled,
               by = c("type", "school", "class")) %>%
    left_join(psu_sampled, by = c("type", "school")) %>%
    select(type, school, class, wt, starts_with("y")) %>%
    mutate(wt = wt / sum(wt) * dplyr::n(),
           across(starts_with("y"), ordered)) %>%
    arrange(type, school, class)

  sampled
}

#' @rdname gen_data_bin_complex
#' @export
gen_data_bin_clust <- gen_data_bin_complex2

gen_data_bin_complex3 <- function(population = make_population(1, seed = NULL),
                                  npsu = 50, n, seed = NULL) {
  # 2-stage stratified cluster sampling
  set.seed(seed)
  if (!missing(n)) npsu <- round(n / (15 + 20 + 25), 0)

  # Weights
  school_info <- population %>%
    select(type, school, class) %>%
    group_by(type, school) %>%
    summarise(
      nstudents = dplyr::n(),
      pr_class_selected = 1 / length(unique(class)),
      .groups = "drop"
    ) %>%
    group_by(type) %>%
    mutate(pr_school_selected = npsu / dplyr::n(),
           prob = pr_school_selected * pr_class_selected,
           wt = 1 / prob)

  # Sampling of PSUs
  psu_sampled <- population %>%
    distinct(type, school) %>%
    group_by(type) %>%
    slice_sample(n = npsu, replace = FALSE)

  # Sampling of classes within PSUs (using SRS)
  classes_sampled <-
    inner_join(population, psu_sampled, by = c("type", "school")) %>%
    distinct(type, school, class) %>%
    group_by(type, school) %>%
    slice_sample(n = 1, replace = FALSE) %>%
    ungroup()

  sampled <-
    inner_join(population, classes_sampled,
               by = c("type", "school", "class")) %>%
    left_join(school_info, by = c("type", "school")) %>%
    select(type, school, class, wt, starts_with("y")) %>%
    mutate(wt = wt / sum(wt) * dplyr::n(),
           across(starts_with("y"), ordered)) %>%
    arrange(type, school, class)

  sampled
}

#' @rdname gen_data_bin_complex
#' @export
gen_data_bin_strcl <- gen_data_bin_complex3
