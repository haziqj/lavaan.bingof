globalVariables(c("i"))

#' Wrapper function for Type 1 and Power simulations
#'
#' @inheritParams txt_mod
#' @inheritParams all_tests
#' @param nsim (integer) The number of simulations to conduct.
#' @param samp_size (integer) The exact sample size for SRS simulations or
#'   stratified sampling; otherwise the *average* sample size for the other
#'   complex sampling methods.
#' @param samp (character) Choose the sampling method for the simulated data.
#'   One of `srs`, `wtd`, `strat`, `clust` or `strcl`.
#' @param simtype (character) Whether this is a `type1` simulation or `power`
#'   simulation.
#' @param starting_seed (integer) The starting random seed.
#' @param ncores (integer) The number of cores to use for parallelisation.
#' @param pop_Sigma (boolean) Should the population value for the multinomial
#'   covariance matrix be used, and not estimated?
#' @param wt (character) Character vector indicating the column name of the
#'   sampling weights to use. Defaults to `NULL`.
#'
#' @return A list of [tibble()]s with the output from [all_tests()].
#' @export
#'
#' @examples
#' \dontrun{
#' library(tidyverse)
#' library(lavaan.bingof)
#' analysis_path <- dirname(rstudioapi::getSourceEditorContext()$path)
#'
#' # Run all scenarios described in manuscript
#' for (sim_type in c("type1", "power")) {
#'   for (samp_method in c("srs", "strat", "clust", "strcl")) {
#'     for (the_samp_size in c(500, 1000, 2000, 3000)) {
#'       for (mod_no in 1:5) {
#'         sim_name <- paste0(samp_method, mod_no, "_n", the_samp_size, "_",
#'                            sim_type)
#'         cat("[", as.character(Sys.time()), "]", "Now running simulation",
#'             sim_name, "\n")
#'         sim <- run_ligof_sims(mod_no, samp_size = the_samp_size,
#'                               samp = samp_method, simtype = sim_type)
#'         invisible(list2env(setNames(list(sim), sim_name), envir = .GlobalEnv))
#'         save(list = sim_name, file = paste0(analysis_path, "/Rsave/",
#'                                             sim_name, ".RData"))
#'       }
#'     }
#'   }
#' }
#' }
run_ligof_sims <- function(model_no = 1, nsim = 1000, samp_size = 1000,
                           samp = c("srs", "wtd", "strat", "clust", "strcl",
                                    "strat2"),
                           simtype = c("type1", "power"), starting_seed = 16423,
                           ncores = parallel::detectCores() - 2,
                           pop_Sigma = FALSE, Sigma2 = NULL, wt = NULL) {

  # Model setup ----------------------------------------------------------------
  mod <- txt_mod(model_no)
  simtype <- match.arg(simtype, c("type1", "power"))
  if (simtype == "type1") H1 <- FALSE
  if (simtype == "power") H1 <- TRUE
  samp <- match.arg(samp, c("srs", "wtd", "strat", "clust", "strcl", "strat2"))

  the_wt <- wt
  if (samp != "srs") {
    # the_wt <- "wt"
    pop <- make_population(model_no, seed = 27324, H1 = H1,
                           Sigma2_attr = isTRUE(pop_Sigma))
    Sigma2pop <- attr(pop, "Sigma2")
  }
  if (isTRUE(pop_Sigma)) {
    Sigma2 <- Sigma2pop
    if (!is.null(Sigma2)) cli::cli_alert_warning("Overriding choice of Sigma2.")
  }

  # Random seeds for replication -----------------------------------------------
  set.seed(starting_seed)
  if (is.null(starting_seed)) {
    the_seeds <- NULL
  } else {
    the_seeds <- 2 ^ (simtype == "power") * model_no * matrix(
      sample(seq_len(100 + nsim ^ 2), size = nsim * 5), ncol = 5
    ) + samp_size
  }

  # Initialise parallel stuff --------------------------------------------------
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  progress <- function(i) setTxtProgressBar(pb, i)
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)
  start_time <- Sys.time()

  # Begin sims -----------------------------------------------------------------
  res <- foreach(
    i = 1:nsim, #.combine = bind_rows,
    .packages = c("dplyr", "tidyr", "purrr", "tibble", "magrittr", "lavaan"),
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
        dat <- gen_data_bin_complex1(population = pop, n = samp_size,
                                     seed = seed_used)
      }
      if (samp == "strat2") {
        # Stratified sampling v2 -----------------------------------------------
        seed_used <- the_seeds[i, 2]
        dat <- gen_data_bin_strat2(model_no, n = samp_size, seed = seed_used,
                                   H1 = H1)
      }
      if (samp == "clust") {
        # Cluster sampling -----------------------------------------------------
        seed_used <- the_seeds[i, 3]
        dat <- gen_data_bin_complex2(population = pop, n = samp_size,
                                     seed = seed_used)
      }
      if (samp == "strcl") {
        # Stratified-cluster sampling ------------------------------------------
        seed_used <- the_seeds[i, 4]
        dat <- gen_data_bin_complex3(population = pop, n = samp_size,
                                     seed = seed_used)
      }
      if (samp == "wtd") {
        # Informative sampling -------------------------------------------------
        seed_used <- the_seeds[i, 5]
        dat <- gen_data_bin_wt(model_no = model_no, n = samp_size,
                               seed = seed_used, H1 = H1)
      }
    }

    fit <- lavaan::sem(model = txt_mod(model_no), data = dat, estimator = "PML",
                       std.lv = TRUE, sampling.weights = the_wt)

    # if (samp == "strat") {
    #   Sigma2 <- create_Sigma2_matrix_complex(fit, strat = dat$type)
    # }
    # if (samp == "clust") {
    #   Sigma2 <- create_Sigma2_matrix_complex(fit, clust = dat$school)
    # }
    # if (samp == "strcl") {
    #   Sigma2 <- create_Sigma2_matrix_complex(fit, strat = dat$type,
    #                                          clust = dat$school)
    # }

    bind_cols(all_tests(fit, sim = i, Sigma2 = Sigma2), seed = seed_used)
  }

  end_time <- Sys.time()
  close(pb)
  stopCluster(cl)

  # Prepare output -------------------------------------------------------------
  class(res) <- "ligof_sims"
  attr(res, "duration") <- difftime(end_time, start_time)
  attr(res, "sim_settings") <- list(
    model_no = model_no,
    nsim = nsim,
    samp_size = samp_size,
    samp = samp,
    simtype = simtype,
    starting_seed = starting_seed,
    ncores = ncores
  )
  res
}

#' Print and summary methods for simulations
#'
#' @name ligof_sims.methods
#'
#' @param x,object The output from [run_ligof_sims()].
#' @param ... Not used.
#'
#' @return The [print()] method displays useful information about the simulation study (number of replications, time taken, model information, etc.). The output of [summary()] is a [tibble()] summarising the rejection rate for an `alpha` level test. For convenience, [summary()] displays the table in the command line interface.
NULL

#' @rdname ligof_sims.methods
#' @export
print.ligof_sims <- function(x, ...) {
  # summary.ligof_sims(x = x, ...)
  print(summary.ligof_sims(x))
}

# print.ligof_sims <- function(x, ...) {
#   list2env(attr(x, "sim_settings"), environment())
#   time_taken <- attr(x, "duration")
#
#   cli::cli_h1("LIGOF simulations")
#
#   cli::cli_text("")
#   cli::cli_text("{.strong Settings}")
#   cli::cli_dl(c("Number of replications" = "{.val {nsim}}"))
#   cli::cli_dl(c("Model" = "{.val {cleanup_model(model_no)}}"))
#   cli::cli_dl(c("Sampling design" = "{.val {cleanup_samp(samp)}}"))
#   cli::cli_dl(c("Sample size" = "{.val {samp_size}}"))
#   cli::cli_text("")
#   cli::cli_text("Simulations completed in {.field {cleanup_duration(time_taken)}}")
# }

cleanup_model <- function(x) {
  paste0(
    x, " (",
    c("1F 5V", "1F 8V", "1F 15V", "2F 10V", "3F 15V")[x],
    ")"
  ) %>% noquote()
}

cleanup_samp <- function(x) {
  res <- "UNKNOWN--NOT CODED YET??"
  if (x == "srs") res <- "Simple random sampling"
  if (x == "strat") res <- "Stratified sampling"
  if (x == "strat2") res <- "Stratified sampling 2"
  if (x == "clust") res <- "Two-stage cluster sampling"
  if (x == "strcl") res <- "Two-stage stratified cluster sampling"
  if (x == "uneqpr") res <- "Two-stage stratified cluster sampling"
  noquote(res)
}

cleanup_duration <- function(x) {
  paste(
    round(as.numeric(x), 1),
    attr(x, "units")
  ) %>% noquote()
}

#' @rdname ligof_sims.methods
#' @param alpha (numeric) The significance level of the test.
#' @export
summary.ligof_sims <- function(object, alpha = 0.05, ...) {
  tab <-
    lapply(object, function(y) if(tibble::is_tibble(y)) { y } else { NULL }) %>%
    bind_rows() %>%
    mutate(alpha_ = pval < alpha,
           name = factor(name, levels = unique(name))) %>%
    group_by(name) %>%
    summarise(n_sims = dplyr::n(),
              n_converged = sum(converged),
              n_rank_def = sum(Omega2_rank < S),
              rej_rate = mean(alpha_[converged & Omega2_rank >= S], na.rm = TRUE),
              mean_X2 = mean(X2[converged & Omega2_rank >= S], na.rm = TRUE),
              mean_df = mean(df[converged & Omega2_rank >= S], na.rm = TRUE),
              .groups = "drop")

  res <- list(
    tab = tab,
    sim_settings = attr(object, "sim_settings"),
    alpha = alpha
  )
  class(res) <- "ligof_sims_summary"
  res
}

#' @export
print.ligof_sims_summary <- function(x, ...) {
  tab <- x$tab
  list2env(x$sim_settings, environment())

  cli::cli_h1("LIGOF simulations summary")
  cli::cli_text("")
  cli::cli_text("Model {.val {cleanup_model(model_no)}} using {.val {tolower(cleanup_samp(samp))}} design (n = {.val {samp_size}})")
  # cli::cli_li("Replications: {.val {nsim}}")
  cli::cli_li("Converged: {.field {tab$n_converged[1]}} / {.val {nsim}}")
  cli::cli_li("Rank deficient: {.field {tab$n_rank_def[1]}} / {.val {nsim}}")
  cli::cli_li("Significance level: {.field {x$alpha}}")

  tab %>%
    select(`Test name` = name, `Rejection rate` = rej_rate,
           `Mean X2 value` = mean_X2, `Mean d.f.` = mean_df) %>%
    kableExtra::kbl(format = "rst", digits = c(1, 3, 2, 2)) %>%
    print()
}
