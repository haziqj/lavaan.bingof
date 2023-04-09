globalVariables(c("i"))

ligof_sims <- function(model.no = 1, nsim = 1000, samp_size = 100,
                       samp = c("srs", "strat", "clust", "strcl"),
                       simtype = c("type1", "power"), starting_seed = 4423,
                       no.cores = parallel::detectCores() - 2) {

  # Model setup
  mod <- txt_mod(model.no)
  simtype <- match.arg(simtype, c("type1", "power"))
  if (simtype == "type1") H1 <- FALSE
  if (simtype == "power") H1 <- TRUE
  samp <- match.arg(samp, c("srs", "strat", "clust", "strcl"))
  the_wt <- NULL
  if (samp != "srs") {
    the_wt <- "wt"
    pop <- make_population(model.no, seed = starting_seed, H1 = H1)
  }

  # Random seeds for replication
  set.seed(starting_seed)
  the_seeds <- model.no * matrix(
    sample(seq_len(100 + nsim ^ 2), size = nsim * 4), ncol = 4
  )

  #
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  progress <- function(i) setTxtProgressBar(pb, i)
  cl <- makeCluster(no.cores)
  registerDoSNOW(cl)

  res <- foreach(
    i = 1:nsim, #.combine = bind_rows,
    .packages = c("tidyverse", "lavaan", "survey"),
    .export = ls(globalenv()),
    .errorhandling = "pass",
    .options.snow = list(progress = progress)
  ) %dopar% {
    if (samp == "srs") {
      # Simple random sampling -------------------------------------------------
      seed_used <- the_seeds[i, 1]
      dat <- gen_data_bin(model.no, n = samp_size, seed = seed_used, H1 = H1)
      svy <- NULL
    } else {
      if (samp == "strat") {
        # Stratified sampling --------------------------------------------------
        seed_used <- the_seeds[i, 2]
        dat <- gen_data_bin_complex1(population = pop, seed = seed_used)
        svy <- svydesign(ids = ~ 1, strata = ~ type, weights = ~ wt, data = dat)
      }
      if (samp == "clust") {
        # Cluster sampling -----------------------------------------------------
        seed_used <- the_seeds[i, 3]
        dat <- gen_data_bin_complex2(population = pop, seed = seed_used)
        svy <- svydesign(ids = ~ school + class, weights = ~ wt, data = dat)
      }
      if (samp == "strcl") {
        # Stratified-cluster sampling ------------------------------------------
        seed_used <- the_seeds[i, 4]
        dat <- gen_data_bin_complex3(population = pop, seed = seed_used)
        svy <- svydesign(ids = ~ school + class, strata = ~ type,
                         weights = ~ wt, data = dat, nest = TRUE)
      }
    }

    fit <- lavaan::sem(model = mod, data = dat, estimator = "PML",
                       std.lv = TRUE, sampling.weights = the_wt)
    bind_cols(all_tests(fit, svy, sim = i), seed = seed_used)
  }

  close(pb)
  stopCluster(cl)

  res
}
