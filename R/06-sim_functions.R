globalVariables(c("i"))

#' Wrapper function for Type 1 and Power simulations
#'
#' @inheritParams txt_mod
#' @param nsim (integer) The number of simulations to conduct.
#' @param samp_size (integer) The exact sample size for SRS simulations or
#'   stratified sampling; otherwise the *average* sample size for the other
#'   complex sampling methods.
#' @param samp (character) Choose the sampling method for the simulated data.
#'   One of `srs`, `strat`, `clust` or `strcl`.
#' @param simtype (character) Whether this is a `type1` simulation or `power`
#'   simulation.
#' @param starting_seed (integer) The starting random seed.
#' @param no.cores (integer) The number of cores to use for parallelisation.
#'
#' @return A list of [tibble()]s with the output from [all_tests()].
#' @export
#'
#' @examples
#' \dontrun{
#' # To run the SRS Type 1 simulations
#' analysis_path <- dirname(rstudioapi::getSourceEditorContext()$path)
#' sim_type <- "type1"
#' for (the_samp_size in c(100, 250, 500, 1000, 2000, 3000)) {
#'   for (mod_no in 1:5) {
#'     sim_name <- paste0("srs", mod_no, "_n", the_samp_size, "_", sim_type)
#'     cat("[", as.character(Sys.time()), "]", "Now running simulation",
#'         sim_name, "\n")
#'     sim <- ligof_sims(mod_no, samp_size = the_samp_size,
#'                       samp = "srs", simtype = sim_type)
#'     list2env(setNames(list(sim), sim_name), envir = .GlobalEnv) %>% invisible
#'     save(list = sim_name, file = paste0(analysis_path, "/Rsave/",
#'                                         sim_name, ".RData"))
#'   }
#' }
#' }
ligof_sims <- function(model_no = 1, nsim = 1000, samp_size = 1000,
                       samp = c("srs", "strat", "clust", "strcl"),
                       simtype = c("type1", "power"), starting_seed = 4423,
                       no.cores = parallel::detectCores() - 2) {

  # Model setup
  mod <- txt_mod(model_no)
  simtype <- match.arg(simtype, c("type1", "power"))
  if (simtype == "type1") H1 <- FALSE
  if (simtype == "power") H1 <- TRUE
  samp <- match.arg(samp, c("srs", "strat", "clust", "strcl"))
  the_wt <- NULL
  if (samp != "srs") {
    the_wt <- "wt"
    pop <- make_population(model_no, seed = starting_seed, H1 = H1)
  }
  if (samp == "strat") npsu <- round(samp_size / 3, 0)
  if (samp == "clust") npsu <- round(samp_size / 21.5, 0)
  if (samp == "strcl") npsu <- round(samp_size / (15 + 20 + 25), 0)

  # Random seeds for replication
  set.seed(starting_seed)
  the_seeds <- 2 ^ (simtype == "power") * model_no * matrix(
    sample(seq_len(100 + nsim ^ 2), size = nsim * 4), ncol = 4
  ) + samp_size

  # Initialise parallel stuff
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  progress <- function(i) setTxtProgressBar(pb, i)
  cl <- makeCluster(no.cores)
  registerDoSNOW(cl)

  start_time <- Sys.time()

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
      dat <- gen_data_bin(model_no, n = samp_size, seed = seed_used, H1 = H1)
      svy <- NULL
    } else {
      if (samp == "strat") {
        # Stratified sampling --------------------------------------------------
        seed_used <- the_seeds[i, 2]
        dat <- gen_data_bin_complex1(population = pop, npsu = npsu,
                                     seed = seed_used)
        svy <- svydesign(ids = ~ 1, strata = ~ type, weights = ~ wt, data = dat)
      }
      if (samp == "clust") {
        # Cluster sampling -----------------------------------------------------
        seed_used <- the_seeds[i, 3]
        dat <- gen_data_bin_complex2(population = pop, npsu = npsu,
                                     seed = seed_used)
        svy <- svydesign(ids = ~ school + class, weights = ~ wt, data = dat)
      }
      if (samp == "strcl") {
        # Stratified-cluster sampling ------------------------------------------
        seed_used <- the_seeds[i, 4]
        dat <- gen_data_bin_complex3(population = pop, npsu = npsu,
                                     seed = seed_used)
        svy <- svydesign(ids = ~ school + class, strata = ~ type,
                         weights = ~ wt, data = dat, nest = TRUE)
      }
    }

    fit <- lavaan::sem(model = mod, data = dat, estimator = "PML",
                       std.lv = TRUE, sampling.weights = the_wt)
    bind_cols(all_tests(fit, svy, sim = i), seed = seed_used)
  }

  end_time <- Sys.time()
  close(pb)
  stopCluster(cl)

  # Prepare output
  class(res) <- "ligof_sims"
  attr(res, "duration") <- difftime(end_time, start_time)
  attr(res, "sim_settings") <- list(
    model_no = model_no,
    nsim = nsim,
    samp_size = samp_size,
    samp = samp,
    simtype = simtype,
    starting_seed = starting_seed,
    no.cores = no.cores
  )
  res
}

#' Print and summary methods for simulations
#'
#' @name ligof_sims.methods
#'
#' @param x,object The output from [ligof_sims()].
#' @param ... Not used.
#'
#' @return To be described.
NULL

#' @rdname ligof_sims.methods
#' @export
print.ligof_sims <- function(x, ...) {
   print(attr(x, "duration"))
}

#' @rdname ligof_sims.methods
#' @export
summary.ligof_sims <- function(object, ...) {
  # cat("Under construction...")
  lapply(object, function(y) if(is_tibble(y)) { y } else { NULL }) %>%
    bind_rows()
}
