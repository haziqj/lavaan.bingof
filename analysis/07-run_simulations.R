library(tidyverse)
library(lavaan.bingof)
analysis_path <- dirname(rstudioapi::getSourceEditorContext()$path)

# All simulations --------------------------------------------------------------
for (sim_type in c("type1", "power")) {
  for (samp_method in c("srs", "strat", "clust", "strcl")) {
    for (the_samp_size in c(500, 1000, 2000, 3000)) {
      for (mod_no in 1:5) {
        sim_name <- paste0(samp_method, mod_no, "_n", the_samp_size, "_", sim_type)
        cat("[", as.character(Sys.time()), "]", "Now running simulation",
            sim_name, "\n")
        sim <- ligof_sims(mod_no, samp_size = the_samp_size, samp = samp_method,
                          simtype = sim_type)
        list2env(setNames(list(sim), sim_name), envir = .GlobalEnv) %>% invisible
        save(list = sim_name, file = paste0(analysis_path, "/Rsave/",
                                            sim_name, ".RData"))
      }
    }
  }
}
