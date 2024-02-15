run_ligof_sims_alt <- function(model_no = 1, nsim = 1000, Sigma2 = NULL,
                           simtype = c("type1", "power"), starting_seed = 16423,
                           ncores = parallel::detectCores() - 2) {

  # Model setup ----------------------------------------------------------------
  mod <- txt_mod(model_no)
  simtype <- match.arg(simtype, c("type1", "power"))
  if (simtype == "type1") H1 <- FALSE
  if (simtype == "power") H1 <- TRUE
  samp <- "strcl2"

  # the_wt <- "wt.norm"
  the_wt <- NULL
  set.seed(123)
  letters <- c(letters, sapply(letters, function(x) paste0(x, letters)))
  pop <-
    tibble(
      type = c("A", "B", "C"),
      fpc1 = c(400, 1000, 600),
      fpc2 = map(fpc1, \(x) { 5 + sample(45, size = x, replace = TRUE)})
    ) |>
    unnest(fpc2) |>
    group_by(type) |>
    mutate(
      sch_id = row_number(),
      tch_id = map(fpc2, \(x) letters[1:x])
    ) |>
    unnest(tch_id) |>
    ungroup()

  pop <- bind_cols(pop, gen_data_bin(model_no = model_no, n = nrow(pop)))

  # Random seeds for replication -----------------------------------------------
  set.seed(starting_seed)

  # Initialise parallel stuff --------------------------------------------------
  pb <- txtProgressBar(min = 0, max = nsim, style = 3)
  progress <- function(i) setTxtProgressBar(pb, i)
  cl <- makeCluster(ncores)
  registerDoSNOW(cl)
  start_time <- Sys.time()

  # Begin sims -----------------------------------------------------------------
  res <- foreach(
    i = 1:nsim, #.combine = bind_rows,
    .packages = c("tidyverse", "lavaan", "lavaan.bingof", "survey"),
    # .export = "pop", #ls(globalenv()),
    .errorhandling = "pass",
    .options.snow = list(progress = progress)
  ) %dopar% {
    dat <- lavaan.bingof:::get_samp(pop = pop)
    fit <- lavaan::sem(model = txt_mod(model_no), data = dat, estimator = "PML",
                       std.lv = TRUE, sampling.weights = the_wt)
    all_tests(fit, sim = i, Sigma2 = Sigma2)
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
    samp_size = 1000,
    samp = samp,
    simtype = simtype,
    starting_seed = starting_seed,
    ncores = ncores
  )
  res
}
