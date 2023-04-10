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
list2env(all_res, envir = .GlobalEnv)

# SRS --------------------------------------------------------------------------

# Type 1 simulations
srs_type1 <- bind_rows(
  # n = 100
  bind_cols(srs1_n100_type1, n = 100, sim = "1F 5V"),
  bind_cols(srs2_n100_type1, n = 100, sim = "1F 8V"),
  bind_cols(srs3_n100_type1, n = 100, sim = "1F 15V"),
  bind_cols(srs4_n100_type1, n = 100, sim = "2F 10V"),
  bind_cols(srs5_n100_type1, n = 100, sim = "3F 15V"),

  # n = 250
  bind_cols(srs1_n250_type1, n = 250, sim = "1F 5V"),
  bind_cols(srs2_n250_type1, n = 250, sim = "1F 8V"),
  bind_cols(srs3_n250_type1, n = 250, sim = "1F 15V"),
  bind_cols(srs4_n250_type1, n = 250, sim = "2F 10V"),
  bind_cols(srs5_n250_type1, n = 250, sim = "3F 15V"),

  # n = 500
  bind_cols(srs1_n500_type1, n = 500, sim = "1F 5V"),
  bind_cols(srs2_n500_type1, n = 500, sim = "1F 8V"),
  bind_cols(srs3_n500_type1, n = 500, sim = "1F 15V"),
  bind_cols(srs4_n500_type1, n = 500, sim = "2F 10V"),
  bind_cols(srs5_n500_type1, n = 500, sim = "3F 15V"),

  # n = 1000
  bind_cols(srs1_n1000_type1, n = 1000, sim = "1F 5V"),
  bind_cols(srs2_n1000_type1, n = 1000, sim = "1F 8V"),
  bind_cols(srs3_n1000_type1, n = 1000, sim = "1F 15V"),
  bind_cols(srs4_n1000_type1, n = 1000, sim = "2F 10V"),
  bind_cols(srs5_n1000_type1, n = 1000, sim = "3F 15V"),

  # n = 2000
  bind_cols(srs1_n2000_type1, n = 2000, sim = "1F 5V"),
  bind_cols(srs2_n2000_type1, n = 2000, sim = "1F 8V"),
  bind_cols(srs3_n2000_type1, n = 2000, sim = "1F 15V"),
  bind_cols(srs4_n2000_type1, n = 2000, sim = "2F 10V"),
  bind_cols(srs5_n2000_type1, n = 2000, sim = "3F 15V"),

  # n = 3000
  bind_cols(srs1_n3000_type1, n = 3000, sim = "1F 5V"),
  bind_cols(srs2_n3000_type1, n = 3000, sim = "1F 8V"),
  bind_cols(srs3_n3000_type1, n = 3000, sim = "1F 15V"),
  bind_cols(srs4_n3000_type1, n = 3000, sim = "2F 10V"),
  bind_cols(srs5_n3000_type1, n = 3000, sim = "3F 15V")
) %>%
  mutate(sim = factor(sim, levels = unique(sim)),
         name = factor(name, levels = rev(unique(name))),
         alpha10 = pval < 0.1,
         alpha5 = pval < 0.05,
         alpha1 = pval < 0.01)

# Power simulations
srs_power <- bind_rows(
  # n = 100
  # bind_cols(srs1_n100_power, n = 100, sim = "1F 5V"),
  # bind_cols(srs2_n100_power, n = 100, sim = "1F 8V"),
  # bind_cols(srs3_n100_power, n = 100, sim = "1F 15V"),
  # bind_cols(srs4_n100_power, n = 100, sim = "2F 10V"),
  # bind_cols(srs5_n100_power, n = 100, sim = "3F 15V"),

  # n = 250
  # bind_cols(srs1_n250_power, n = 250, sim = "1F 5V"),
  # bind_cols(srs2_n250_power, n = 250, sim = "1F 8V"),
  # bind_cols(srs3_n250_power, n = 250, sim = "1F 15V"),
  # bind_cols(srs4_n250_power, n = 250, sim = "2F 10V"),
  # bind_cols(srs5_n250_power, n = 250, sim = "3F 15V"),

  # n = 500
  bind_cols(srs1_n500_power, n = 500, sim = "1F 5V"),
  # bind_cols(srs2_n500_power, n = 500, sim = "1F 8V"),
  # bind_cols(srs3_n500_power, n = 500, sim = "1F 15V"),
  # bind_cols(srs4_n500_power, n = 500, sim = "2F 10V"),
  # bind_cols(srs5_n500_power, n = 500, sim = "3F 15V"),

  # n = 1000
  bind_cols(srs1_n1000_power, n = 1000, sim = "1F 5V"),
  # bind_cols(srs2_n1000_power, n = 1000, sim = "1F 8V"),
  # bind_cols(srs3_n1000_power, n = 1000, sim = "1F 15V"),
  # bind_cols(srs4_n1000_power, n = 1000, sim = "2F 10V"),
  # bind_cols(srs5_n1000_power, n = 1000, sim = "3F 15V"),

  # n = 2000
  bind_cols(srs1_n2000_power, n = 2000, sim = "1F 5V"),
  # bind_cols(srs2_n2000_power, n = 2000, sim = "1F 8V"),
  # bind_cols(srs3_n2000_power, n = 2000, sim = "1F 15V"),
  # bind_cols(srs4_n2000_power, n = 2000, sim = "2F 10V"),
  # bind_cols(srs5_n2000_power, n = 2000, sim = "3F 15V"),

  # n = 3000
  bind_cols(srs1_n3000_power, n = 3000, sim = "1F 5V"),
  # bind_cols(srs2_n3000_power, n = 3000, sim = "1F 8V"),
  # bind_cols(srs3_n3000_power, n = 3000, sim = "1F 15V"),
  # bind_cols(srs4_n3000_power, n = 3000, sim = "2F 10V"),
  # bind_cols(srs5_n3000_power, n = 3000, sim = "3F 15V")
) %>%
  mutate(sim = factor(sim, levels = unique(sim)),
         name = factor(name, levels = rev(unique(name))),
         alpha10 = pval < 0.1,
         alpha5 = pval < 0.05,
         alpha1 = pval < 0.01)

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
    # geom_ribbon(aes(ymin = rej_rate10 - 1.96 * se10,
    #                 ymax = rej_rate10 + 1.96 * se10)) +
    facet_wrap(. ~ sim, ncol = 3) +
    scale_x_continuous(breaks = unique(srs_power$n)) +
    scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 11)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = c(0.85, 0.24)) +
    labs(x = "Sample size (n)", y = "Rejection proportion", col = NULL,
         shape = NULL, title = as.expression(bquote(
           .(plot_title)~"("*alpha~"="~.(iprior::dec_plac(alpha/100, 2))*")"
         ))) +
    scale_colour_viridis_d(option = "turbo")
}

# SRS results ------------------------------------------------------------------

# Type I error results
srs_type1_res <-
  srs_type1 %>%
  group_by(name, sim, n) %>%
  summarise(n_sims = n(),
            n_converged = sum(converged),
            max_rank = max(Omega2_rank),
            n_rank_def = sum(Omega2_rank < max_rank),
            rej_rate10 = mean(alpha10[converged], na.rm = TRUE),
            # se10  = sd(alpha10, na.rm = TRUE) / n(),
            rej_rate5 = mean(alpha5[converged], na.rm = TRUE),
            # se5  = sd(alpha5, na.rm = TRUE) / n(),
            rej_rate1 = mean(alpha1[converged], na.rm = TRUE),
            # se1  = sd(alpha1, na.rm = TRUE) / n(),
            .groups = "drop")

# Power results
srs_power_res <-
  srs_power %>%
  group_by(name, sim, n) %>%
  summarise(n_sims = n(),
            n_converged = sum(converged),
            max_rank = max(Omega2_rank),
            n_rank_def = sum(Omega2_rank < max_rank),
            rej_rate10 = mean(alpha10[converged], na.rm = TRUE),
            # se10  = sd(alpha10, na.rm = TRUE) / n(),
            rej_rate5 = mean(alpha5[converged], na.rm = TRUE),
            # se5  = sd(alpha5, na.rm = TRUE) / n(),
            rej_rate1 = mean(alpha1[converged], na.rm = TRUE),
            # se1  = sd(alpha1, na.rm = TRUE) / n(),
            .groups = "drop")

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
     p_srs_e, p_srs_f, file = "simres_srs.RData")

# Complex sampling -------------------------------------------------------------
strat_type1 <- bind_rows(
  bind_cols(strat1_type1, sim = "1F 5V"),
  bind_cols(strat2_type1, sim = "1F 8V"),
  bind_cols(strat3_type1, sim = "1F 15V"),
  bind_cols(strat4_type1, sim = "2F 10V"),
  bind_cols(strat5_type1, sim = "3F 15V")
) %>%
  mutate(sim = factor(sim, levels = unique(sim)),
         name = factor(name, levels = rev(unique(name))),
         alpha10 = pval < 0.1,
         alpha5 = pval < 0.05,
         alpha1 = pval < 0.01)

strat_power <- bind_rows(
  bind_cols(strat1_power, sim = "1F 5V"),
  bind_cols(strat2_power, sim = "1F 8V"),
  bind_cols(strat3_power, sim = "1F 15V"),
  bind_cols(strat4_power, sim = "2F 10V"),
  bind_cols(strat5_power, sim = "3F 15V")
) %>%
  mutate(sim = factor(sim, levels = unique(sim)),
         name = factor(name, levels = rev(unique(name))),
         alpha10 = pval < 0.1,
         alpha5 = pval < 0.05,
         alpha1 = pval < 0.01)

clust_type1 <- bind_rows(
  bind_cols(clust1_type1, sim = "1F 5V"),
  bind_cols(clust2_type1, sim = "1F 8V"),
  bind_cols(clust3_type1, sim = "1F 15V"),
  bind_cols(clust4_type1, sim = "2F 10V"),
  bind_cols(clust5_type1, sim = "3F 15V")
) %>%
  mutate(sim = factor(sim, levels = unique(sim)),
         name = factor(name, levels = rev(unique(name))),
         alpha10 = pval < 0.1,
         alpha5 = pval < 0.05,
         alpha1 = pval < 0.01)

clust_power <- bind_rows(
  bind_cols(clust1_power, sim = "1F 5V"),
  bind_cols(clust2_power, sim = "1F 8V"),
  bind_cols(clust3_power, sim = "1F 15V"),
  bind_cols(clust4_power, sim = "2F 10V"),
  bind_cols(clust5_power, sim = "3F 15V")
) %>%
  mutate(sim = factor(sim, levels = unique(sim)),
         name = factor(name, levels = rev(unique(name))),
         alpha10 = pval < 0.1,
         alpha5 = pval < 0.05,
         alpha1 = pval < 0.01)

strcl_type1 <- bind_rows(
  bind_cols(strcl1_type1, sim = "1F 5V"),
  bind_cols(strcl2_type1, sim = "1F 8V"),
  bind_cols(strcl3_type1, sim = "1F 15V"),
  bind_cols(strcl4_type1, sim = "2F 10V"),
  bind_cols(strcl5_type1, sim = "3F 15V")
) %>%
  mutate(sim = factor(sim, levels = unique(sim)),
         name = factor(name, levels = rev(unique(name))),
         alpha10 = pval < 0.1,
         alpha5 = pval < 0.05,
         alpha1 = pval < 0.01)

strcl_power <- bind_rows(
  bind_cols(strcl1_power, sim = "1F 5V"),
  bind_cols(strcl2_power, sim = "1F 8V"),
  bind_cols(strcl3_power, sim = "1F 15V"),
  bind_cols(strcl4_power, sim = "2F 10V"),
  bind_cols(strcl5_power, sim = "3F 15V")
) %>%
  mutate(sim = factor(sim, levels = unique(sim)),
         name = factor(name, levels = rev(unique(name))),
         alpha10 = pval < 0.1,
         alpha5 = pval < 0.05,
         alpha1 = pval < 0.01)

complex_type1 <- bind_rows(
  bind_cols(strat_type1, sampling = "Stratified"),
  bind_cols(clust_type1, sampling = "Cluster"),
  bind_cols(strcl_type1, sampling = "Strat-clust")
) %>%
  mutate(sampling = factor(sampling, levels = c("Stratified", "Cluster", "Strat-clust")))

complex_power <- bind_rows(
  bind_cols(strat_power, sampling = "Stratified"),
  bind_cols(clust_power, sampling = "Cluster"),
  bind_cols(strcl_power, sampling = "Strat-clust")
) %>%
  mutate(sampling = factor(sampling, levels = c("Stratified", "Cluster", "Strat-clust")))

complex_plot <- function(x = complex_type1_res, alpha = 10, dashed_line = TRUE,
                         plot_title = "Type I errors") {
  var_name <- paste0("rej_rate", alpha)

  p <-
    ggplot(x, aes(sampling, .data[[var_name]], fill = name,
                  alpha = n_rank_def / 10)) +
    geom_bar(stat = "identity", position = "dodge")

  if (isTRUE(dashed_line)) {
    p <- p +
      geom_hline(yintercept = alpha / 100, linetype = "dashed", col = "grey50")

  }
  p +
    facet_wrap(. ~ sim, ncol = 3) +
    scale_alpha("% rank def.", range = c(1, 0.3)) +
    theme(legend.position = c(0.85, 0.24)) +
    labs(x = "Complex sampling method", y = "Rejection proportion", fill = NULL,
         shape = NULL, title = as.expression(bquote(
           .(plot_title)~"("*alpha~"="~.(iprior::dec_plac(alpha/100, 2))*")"
         ))) +
    guides(fill = guide_legend(ncol = 2), alpha = guide_legend(ncol = 2)) +
    scale_fill_viridis_d(option = "turbo")
}

# Complex sampling results -----------------------------------------------------
complex_type1_res <-
  complex_type1 %>%
  group_by(name, sim, sampling) %>%
  summarise(n_sims = n(),
            n_converged = sum(converged),
            max_rank = max(Omega2_rank),
            n_rank_def = sum(Omega2_rank < max_rank),
            rej_rate10 = mean(alpha10[converged], na.rm = TRUE),
            # se10  = sd(alpha10, na.rm = TRUE) / n(),
            rej_rate5 = mean(alpha5[converged], na.rm = TRUE),
            # se5  = sd(alpha5, na.rm = TRUE) / n(),
            rej_rate1 = mean(alpha1[converged], na.rm = TRUE),
            # se1  = sd(alpha1, na.rm = TRUE) / n(),
            .groups = "drop")

# Power results
complex_power_res <-
  complex_power %>%
  group_by(name, sim, sampling) %>%
  summarise(n_sims = n(),
            n_converged = sum(converged),
            max_rank = max(Omega2_rank),
            n_rank_def = sum(Omega2_rank < max_rank),
            rej_rate10 = mean(alpha10[converged], na.rm = TRUE),
            # se10  = sd(alpha10, na.rm = TRUE) / n(),
            rej_rate5 = mean(alpha5[converged], na.rm = TRUE),
            # se5  = sd(alpha5, na.rm = TRUE) / n(),
            rej_rate1 = mean(alpha1[converged], na.rm = TRUE),
            # se1  = sd(alpha1, na.rm = TRUE) / n(),
            .groups = "drop")

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
