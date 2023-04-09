library(tidyverse)
library(lavaan.bingof)
analysis_path <- dirname(rstudioapi::getSourceEditorContext()$path)

# SRS simulations --------------------------------------------------------------

# Type 1 errors
sim_type <- "type1"
for (the_samp_size in c(100, 250, 500, 1000, 2000, 3000)) {
  for (mod_no in 1:5) {
    sim_name <- paste0("srs", mod_no, "_n", the_samp_size, "_", sim_type)
    cat("[", as.character(Sys.time()), "]", "Now running simulation",
        sim_name, "\n")
    sim <- lavaan.bingof:::ligof_sims(mod_no, samp_size = the_samp_size,
                                      samp = "srs", simtype = sim_type)
    list2env(setNames(list(sim), sim_name), envir = .GlobalEnv) %>% invisible
    save(list = sim_name, file = paste0(analysis_path, "/Rsave/",
                                        sim_name, ".RData"))
  }
}

# Power
sim_type <- "power"
for (the_samp_size in c(100, 250, 500, 1000, 2000, 3000)) {
  for (mod_no in 1:5) {
    sim_name <- paste0("srs", mod_no, "_n", the_samp_size, "_", sim_type)
    cat("[", as.character(Sys.time()), "]", "Now running simulation",
        sim_name, "\n")
    sim <- lavaan.bingof:::ligof_sims(mod_no, samp_size = the_samp_size,
                                      samp = "srs", simtype = sim_type)
    list2env(setNames(list(sim), sim_name), envir = .GlobalEnv) %>% invisible
    save(list = sim_name, file = paste0(analysis_path, "/Rsave/",
                                        sim_name, ".RData"))
  }
}

# Complex sampling simulations -------------------------------------------------

# Type 1 errors
sim_type <- "type1"
for (samp_method in c("strat", "clust", "strcl")) {
  for (mod_no in 1:5) {
    sim_name <- paste0(samp_method, mod_no, "_", sim_type)
    cat("[", as.character(Sys.time()), "]", "Now running simulation",
        sim_name, "\n")
    sim <- lavaan.bingof:::ligof_sims(mod_no, samp = samp_method,
                                      simtype = sim_type)
    list2env(setNames(list(sim), sim_name), envir = .GlobalEnv) %>% invisible
    save(list = sim_name, file = paste0(analysis_path, "/Rsave/", sim_name,
                                        ".RData"))
  }
}

# Power
sim_type <- "power"
for (samp_method in c("strat", "clust", "strcl")) {
  for (mod_no in 1:5) {
    sim_name <- paste0(samp_method, mod_no, "_", sim_type)
    cat("[", as.character(Sys.time()), "]", "Now running simulation",
        sim_name, "\n")
    sim <- lavaan.bingof:::ligof_sims(mod_no, samp = samp_method,
                                      simtype = sim_type)
    list2env(setNames(list(sim), sim_name), envir = .GlobalEnv) %>% invisible
    save(list = sim_name, file = paste0(analysis_path, "/Rsave/", sim_name,
                                        ".RData"))
  }
}
