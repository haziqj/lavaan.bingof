library(tidyverse)
theme_set(theme_bw())

# Load saved simulations -------------------------------------------------------
rm(list = ls())
sim_saved <- paste0("analysis/Rsave/", dir("analysis/Rsave/"))
for (i in sim_saved) load(i)
rm(list = c("sim_saved", "i"))
all_res <- mget(ls())

# Combine function -------------------------------------------------------------
my_comb_fn <- function(x) {
  lapply(x, function(y) if(is_tibble(y)) { y } else { NULL }) %>%
    bind_rows()
}
all_res <- lapply(all_res, my_comb_fn)  # clean all results

grab_sims <- function(x = all_res, samp = "srs", type = "type1") {
  mod_names <- c("1F 5V", "1F 8V", "1F 15V", "2F 10V", "3F 15V")

  type_of_sampling <- grepl(samp, names(x))
  type_of_analysis <- grepl(type, names(x))
  res <- list(NULL)

  for (n in c(500, 1000, 2000, 3000)) {
    type_of_n <- grepl(n, names(x))
    ind <- which(type_of_sampling & type_of_analysis & type_of_n)
    tmp <- NULL
    for (i in seq_along(ind)) {
      tmp <- bind_rows(tmp, bind_cols(x[[ind[i]]], n = n, sim = mod_names[i]))
    }
    res <- c(res, list(tmp))
  }

  do.call("bind_rows", res) %>%
    mutate(sim = factor(sim, levels = unique(sim)),
           name = factor(name, levels = unique(name)),
           alpha10 = pval < 0.1,
           alpha5 = pval < 0.05,
           alpha1 = pval < 0.01)
}

summarise_sims <- function(samp = "srs", type = "type1") {
  grab_sims(samp = samp, type = type) %>%
    group_by(name, sim, n) %>%
    summarise(n_sims = n(),
              n_converged = sum(converged),
              n_rank_def = sum(Omega2_rank < S),
              rej_rate10 = mean(alpha10[converged], na.rm = TRUE),
              rej_rate5 = mean(alpha5[converged], na.rm = TRUE),
              rej_rate1 = mean(alpha1[converged], na.rm = TRUE),
              .groups = "drop")
}

# Plot functions ---------------------------------------------------------------
srs_plot <- function(x = srs_type1_res, alpha = 10, dashed_line = TRUE,
                     plot_title = "Type I errors") {
  var_name <- paste0("rej_rate", alpha)

  p <- ggplot(x, aes(n, .data[[var_name]], col = name, shape = name))
  if (isTRUE(dashed_line)) {
    p <- p + geom_hline(yintercept = alpha / 100, linetype = "dashed",
                        col = "grey50")
  }
  p +
    geom_point() +
    geom_line() +
    facet_wrap(. ~ sim, ncol = 3) +
    scale_x_continuous(breaks = unique(srs_power_res$n)) +
    scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 11)) +
    scale_alpha("% rank def.", range = c(1, 0.3)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = c(0.85, 0.2)) +
    labs(x = "Sample size (n)", y = "Rejection proportion", col = NULL,
         shape = NULL, title = as.expression(bquote(
           .(plot_title)~"("*alpha~"="~.(iprior::dec_plac(alpha/100, 2))*")"
         ))) +
    guides(col = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
    scale_colour_viridis_d(option = "turbo", direction = -1)
}

# SRS results ------------------------------------------------------------------
srs_type1_res <- summarise_sims("srs", "type1")
srs_power_res <- summarise_sims("srs", "power")

p_srs_a <- srs_plot(srs_type1_res, alpha = 10)
p_srs_b <- srs_plot(srs_type1_res, alpha = 5)
p_srs_c <- srs_plot(srs_type1_res, alpha = 1)
p_srs_d <- srs_plot(srs_power_res, alpha = 10, dashed_line = FALSE,
                    plot_title = "Power")
p_srs_e <- srs_plot(srs_power_res, alpha = 5, dashed_line = FALSE,
                    plot_title = "Power")
p_srs_f <- srs_plot(srs_power_res, alpha = 1, dashed_line = FALSE,
                    plot_title = "Power")

save(srs_type1_res, srs_power_res, p_srs_a, p_srs_b, p_srs_c, p_srs_d,
     p_srs_e, p_srs_f, file = "analysis/simres_srs.RData")

# Complex sampling results -----------------------------------------------------
complex_type1_res <- bind_rows(
  bind_cols(summarise_sims("strat", "type1"), sampling = "Stratified"),
  bind_cols(summarise_sims("clust", "type1"), sampling = "Cluster"),
  bind_cols(summarise_sims("strcl", "type1"), sampling = "Strat-clust")
) %>%
  mutate(sampling = factor(sampling, levels = c("Stratified", "Cluster",
                                                "Strat-clust")))

complex_power_res <- bind_rows(
  bind_cols(summarise_sims("strat", "power"), sampling = "Stratified"),
  bind_cols(summarise_sims("clust", "power"), sampling = "Cluster"),
  bind_cols(summarise_sims("strcl", "power"), sampling = "Strat-clust")
) %>%
  mutate(sampling = factor(sampling, levels = c("Stratified", "Cluster",
                                                "Strat-clust")))

complex_plot <- function(x = complex_type1_res, alpha = 10, dashed_line = TRUE,
                         plot_title = "Type I errors") {
  var_name <- paste0("rej_rate", alpha)

  p <-
    x %>%
    mutate(n = factor(n, labels = paste0("n =\n", unique(x$n)))) %>%
    ggplot(aes(n, .data[[var_name]], fill = name, alpha = n_rank_def / 10)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.9)

  if (isTRUE(dashed_line)) {
    p <- p +
      geom_hline(aes(yintercept = alpha / 100, linetype = "Nominal\nrej. level"),
                 col = "grey50") +
      scale_linetype_manual(NULL, values = "dashed")
  }
  p +
    facet_wrap(. ~ sim, ncol = 3) +
    scale_alpha("% rank\ndef.", range = c(1, 0.3)) +
    theme(legend.position = "bottom") +
    labs(x = "Sample size", y = "Rejection proportion", fill = NULL,
         shape = NULL, title = as.expression(bquote(
           .(plot_title)~"("*alpha~"="~.(iprior::dec_plac(alpha/100, 2))*")"
         ))) +
    guides(fill = guide_legend(ncol = 3, order = 1),
           alpha = guide_legend(ncol = 2, order = 2)) +
    scale_fill_viridis_d(option = "turbo", direction = -1) +
    facet_grid(sim ~ sampling)
}

p_complex_a <- complex_plot(complex_type1_res, alpha = 10)
p_complex_b <- complex_plot(complex_type1_res, alpha = 5)
p_complex_c <- complex_plot(complex_type1_res, alpha = 1)
p_complex_d <- complex_plot(complex_power_res, alpha = 10, dashed_line = FALSE,
                            plot_title = "Power")
p_complex_e <- complex_plot(complex_power_res, alpha = 5, dashed_line = FALSE,
                            plot_title = "Power")
p_complex_f <- complex_plot(complex_power_res, alpha = 1, dashed_line = FALSE,
                            plot_title = "Power")

save(complex_type1_res, complex_power_res, p_complex_a, p_complex_b,
     p_complex_c, p_complex_d, p_complex_e, p_complex_f,
     file = "analysis/simres_complex.RData")

files_to_copy <- c("simres_srs.RData", "simres_complex.RData")
file.copy(from = paste0("analysis/", files_to_copy),
          to = paste0("vignettes/articles/", files_to_copy), overwrite = TRUE)
rm(list = ls())
