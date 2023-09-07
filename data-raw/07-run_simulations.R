library(tidyverse)
library(lavaan.bingof)
library(survey)
analysis_path <- dirname(rstudioapi::getSourceEditorContext()$path)

# All simulations --------------------------------------------------------------
for (sim_type in c("type1")) {
  for (samp_method in c("strat")) {
    for (the_samp_size in c(2500)) {
      for (mod_no in 1:5) {
        sim_name <- paste0(samp_method, mod_no, "_n", the_samp_size, "_",
                           sim_type)
        cat("[", as.character(Sys.time()), "]", "Now running simulation",
            sim_name, "\n")
        sim <- run_ligof_sims(mod_no, samp_size = the_samp_size, nsim = 250,
                              samp = samp_method, simtype = sim_type,
                              pop_Sigma = FALSE, bootstrap = FALSE)
        invisible(list2env(setNames(list(sim), sim_name), envir = .GlobalEnv))
        # save(list = sim_name, file = paste0(analysis_path, "/Rsave/",
        #                                     sim_name, ".RData"))
      }
    }
  }
}

# Create data file -------------------------------------------------------------

# Uncomment these lines if loading first is necessary
rm(list = ls())
analysis_path <- dirname(rstudioapi::getSourceEditorContext()$path)
sim_saved <-
  paste0(analysis_path, "/Rsave/") %>%
  paste0(., dir(.))
for (i in sim_saved) load(i)
rm(list = c("sim_saved", "i", "analysis_path"))

# Clean results
all_res <-
  mget(ls()) %>%
  lapply(., FUN = \(x) {
    lapply(x, \(y) if(is_tibble(y)) { y } else { NULL }) %>%
      bind_rows()
  })

usethis::use_data(all_res, overwrite = TRUE, compress = "xz")
