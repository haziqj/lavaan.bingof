library(tidyverse)
theme_set(theme_bw())
library(jcolors)
library(lavaan.bingof)  # needs all_res data. be sure to build pkg
# data(all_res, package = "lavaan.bingof")

grab_sims <- function(x = all_res, samp = "srs", type = "type1",
                      the_n = c(500, 1000, 2000, 3000, 5000, 10000)) {
  mod_names <- c("1F 5V", "1F 8V", "1F 15V", "2F 10V", "3F 15V")

  type_of_sampling <- grepl(samp, names(x))
  type_of_analysis <- grepl(type, names(x))
  res <- list(NULL)

  for (n in the_n) {
    type_of_n <- grepl(paste0(n, "_"), names(x))
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
              crit10 = qnorm(0.975) * sqrt(rej_rate10 * (1 - rej_rate10) /
                                             n_converged),
              rej_rate5 = mean(alpha5[converged], na.rm = TRUE),
              crit5 = qnorm(0.975) * sqrt(rej_rate5 * (1 - rej_rate5) /
                                            n_converged),
              rej_rate1 = mean(alpha1[converged], na.rm = TRUE),
              crit1 = qnorm(0.975) * sqrt(rej_rate1 * (1 - rej_rate1) /
                                            n_converged),
              .groups = "drop")
}

# SRS results ------------------------------------------------------------------
res_srs_type1 <- summarise_sims("srs", "type1")
res_srs_power <- summarise_sims("srs", "power")

srs_plot <- function(x = res_srs_type1, alpha = 10, dashed_line = TRUE,
                     plot_title = "Type I errors",
                     exclude_tests = c("RSS,MM3", "Multn,MM3")) {
  var_name <- paste0("rej_rate", alpha)
  crit_name <- paste0("crit", alpha)

  x <- x %>%
    filter(!name %in% exclude_tests)

  p <- ggplot(x, aes(n, .data[[var_name]], col = name, shape = name))
  if (isTRUE(dashed_line)) {
    p <- p + geom_hline(yintercept = alpha / 100, linetype = "dashed",
                        col = "grey50")
  }
  p +
    geom_ribbon(aes(ymin = .data[[var_name]] - .data[[crit_name]],
                    ymax = .data[[var_name]] + .data[[crit_name]],
                    fill = name), col = NA, alpha = 0.1) +
    geom_point() +
    geom_line() +
    facet_wrap(. ~ sim, ncol = 3) +
    scale_x_continuous(breaks = unique(res_srs_type1$n)) +
    scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 11, 16, 17, 15, 3, 7)) +
    scale_alpha("% rank def.", range = c(1, 0.3)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = c(0.85, 0.2)) +
    labs(x = "Sample size (n)", y = "Rejection proportion", col = NULL,
         fill = NULL,
         shape = NULL, title = as.expression(bquote(
           .(plot_title)~"("*alpha~"="~.(iprior::dec_plac(alpha/100, 2))*")"
         ))) +
    guides(col = guide_legend(ncol = 1), shape = guide_legend(ncol = 1)) +
    # scale_colour_viridis_d(option = "turbo", direction = -1) +
    # scale_fill_viridis_d(option = "turbo", direction = -1) +
    jcolors::scale_colour_jcolors(palette = "pal8") +
    jcolors::scale_fill_jcolors(palette = "pal8")
}

srs_plot2 <- function(x = res_srs_type1, alpha = 5, dashed_line = TRUE,
                      plot_title = "Type I errors",
                      exclude_tests = c("RSS,MM3", "Multn,MM3")) {
  var_name <- paste0("rej_rate", alpha)
  crit_name <- paste0("crit", alpha)

  x <- x %>%
    filter(!name %in% exclude_tests) %>%
    mutate(sim = factor(sim, labels = c("1*F~5*V", "1*F~8*V", "1*F~15*V",
                                        "2*F~10*V", "3*F~15*V")),
           name = factor(name, levels = rev(levels(name))),
           ok = !(!is.na(.data[[var_name]]) & !(abs(.data[[var_name]] - alpha / 100) < .data[[crit_name]])),
           N = factor(n))

  p <- ggplot(x, aes(.data[[var_name]], name, col = N, shape = N, alpha = ok)) +
    geom_vline(aes(xintercept = alpha / 100), linetype = "dashed") +
    geom_pointrange(aes(xmin = .data[[var_name]] - .data[[crit_name]],
                        xmax = .data[[var_name]] + .data[[crit_name]]),
                    position = position_dodge(width = 0.5)) +
    facet_grid(. ~ sim, labeller = label_parsed) +
    jcolors::scale_colour_jcolors(palette = "pal8") +
    scale_alpha_manual(values = c(0.4, 1)) +
    # scale_x_continuous(breaks = alpha / 100) +
    scale_shape_manual(values = c(16, 17, 15, 4, 1, 2)) +
    labs(x = "Rejection rate", y = NULL, alpha = "Within\n95% interval",
         shape = "Sample size", col = "Sample size") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  p
}

p_srs_a <- srs_plot2(res_srs_type1, alpha = 10)
p_srs_b <- srs_plot2(res_srs_type1, alpha = 5)
p_srs_c <- srs_plot2(res_srs_type1, alpha = 1)
p_srs_d <- srs_plot(res_srs_power, alpha = 10, dashed_line = FALSE,
                    plot_title = "Power") +
  scale_colour_viridis_d(option = "turbo", direction = -1) +
  scale_fill_viridis_d(option = "turbo", direction = -1)
p_srs_e <- srs_plot(res_srs_power, alpha = 5, dashed_line = FALSE,
                    plot_title = "Power") +
  scale_colour_viridis_d(option = "turbo", direction = -1) +
  scale_fill_viridis_d(option = "turbo", direction = -1)
p_srs_f <- srs_plot(res_srs_power, alpha = 1, dashed_line = FALSE,
                    plot_title = "Power") +
  scale_colour_viridis_d(option = "turbo", direction = -1) +
  scale_fill_viridis_d(option = "turbo", direction = -1)

usethis::use_data(res_srs_type1, overwrite = TRUE)
usethis::use_data(res_srs_power, overwrite = TRUE)

# Complex sampling results -----------------------------------------------------
res_complex_type1 <- bind_rows(
  bind_cols(summarise_sims("strat", "type1"), sampling = "Stratified"),
  bind_cols(summarise_sims("clust", "type1"), sampling = "Cluster"),
  bind_cols(summarise_sims("strcl", "type1"), sampling = "Strat-clust")
) %>%
  mutate(sampling = factor(sampling, levels = c("Stratified", "Cluster",
                                                "Strat-clust")))

res_complex_power <- bind_rows(
  bind_cols(summarise_sims("strat", "power"), sampling = "Stratified"),
  bind_cols(summarise_sims("clust", "power"), sampling = "Cluster"),
  bind_cols(summarise_sims("strcl", "power"), sampling = "Strat-clust")
) %>%
  mutate(sampling = factor(sampling, levels = c("Stratified", "Cluster",
                                                "Strat-clust")))

complex_plot <- function(x = res_complex_type1, alpha = 10, dashed_line = TRUE,
                         plot_title = "Type I errors",
                         exclude_tests = c("RSS,MM3", "Multn,MM3"),
                         exclude_sims = NULL) {
  var_name <- paste0("rej_rate", alpha)
  crit_name <- paste0("crit", alpha)
  nsim <- res_complex_type1$n_sims[1]
  x <- x %>%
    filter(!name %in% exclude_tests) %>%
    filter(!sim %in% exclude_sims)

  p <-
    x %>%
    mutate(n = factor(n, labels = paste0("n =\n", unique(x$n)))) %>%
    ggplot(aes(n, .data[[var_name]], fill = name,
               alpha = n_rank_def / nsim * 100
               )) +
    geom_bar(stat = "identity", position = "dodge", width = 0.9) +
    geom_errorbar(aes(ymin = .data[[var_name]] - .data[[crit_name]],
                      ymax = .data[[var_name]] + .data[[crit_name]]),
                  col = "black", width = 0.2, alpha = 1,
                  position = position_dodge(width = 0.9))

  if (isTRUE(dashed_line)) {
    p <- p +
      geom_hline(aes(yintercept = alpha / 100, linetype = "Nominal\nrej. level"),
                 col = "grey50") +
      scale_linetype_manual(NULL, values = "dashed")
  }
  p +
    facet_wrap(. ~ sim, ncol = 3) +
    scale_alpha("% rank\ndef.", range = c(1, 0.5)) +
    theme(legend.position = "bottom") +
    labs(x = "Sample size", y = "Rejection proportion", fill = NULL,
         shape = NULL, title = as.expression(bquote(
           .(plot_title)~"("*alpha~"="~.(iprior::dec_plac(alpha/100, 2))*")"
         ))) +
    guides(fill = guide_legend(nrow = 2, order = 1),
           alpha = guide_legend(ncol = 3, order = 2)) +
    scale_fill_viridis_d(option = "turbo", direction = -1) +
    scale_colour_viridis_d(option = "turbo", direction = -1) +
    # scale_fill_jcolors() +
    facet_grid(sim ~ sampling)
}

p_complex_a <- complex_plot(res_complex_type1, alpha = 10) +
  coord_cartesian(ylim = c(0, 0.2))
p_complex_b <- complex_plot(res_complex_type1, alpha = 5) +
  coord_cartesian(ylim = c(0, 0.2))
p_complex_c <- complex_plot(res_complex_type1, alpha = 1) +
  coord_cartesian(ylim = c(0, 0.05))
p_complex_d <- complex_plot(res_complex_power, alpha = 10, dashed_line = FALSE,
                            plot_title = "Power")
p_complex_e <- complex_plot(res_complex_power, alpha = 5, dashed_line = FALSE,
                            plot_title = "Power")
p_complex_f <- complex_plot(res_complex_power, alpha = 1, dashed_line = FALSE,
                            plot_title = "Power")

usethis::use_data(res_complex_type1, overwrite = TRUE)
usethis::use_data(res_complex_power, overwrite = TRUE)

# Distribution of test statistics ----------------------------------------------
plot_X2_dens <- function(samp = "strcl") {
  dat <- grab_sims(samp = samp)
  plot_df2 <-
    dat %>%
    group_by(name, sim) %>%
    summarise(mdf = mean(df),
              sdf = sd(df),
              minx = min(X2),
              maxx = max(X2),
              .groups = "drop") %>%
    rename(df = mdf) %>%
    rowwise() %>%
    mutate(xx = list(seq(minx, maxx, length = 200)),
           yy = list(dchisq(xx, df = df))) %>%
    unnest_longer(col = c(xx, yy))


  dat %>%
    mutate(n = paste0("n = ", n),
           n = factor(n, levels = paste0("n = ", c(500, 1000, 2000, 3000, 5000,
                                                   10000)))) %>%
    ggplot(aes(X2, y = stat(density), group = n, fill = n)) +
    geom_histogram(position = "identity", alpha = 0.5, col = NA, bins = 30) +
    # geom_density(alpha = 0.5) +
    # geom_freqpoly() +
    geom_line(data = plot_df2, aes(xx, yy), col = "black", inherit.aes = FALSE) +
    facet_wrap(name ~ sim, scales = "free", ncol = 5) +
    jcolors::scale_colour_jcolors(palette = "pal7", name = NULL) +
    jcolors::scale_fill_jcolors(palette = "pal7", name = NULL) +
    theme(legend.position = "bottom") +
    guides(fill = guide_legend(nrow = 1))
}

p_hist_a <- plot_X2_dens("srs")
p_hist_b <- plot_X2_dens("strat")
p_hist_c <- plot_X2_dens("clust")
p_hist_d <- plot_X2_dens("strcl")

# dat %>%
#   ggplot(aes(x = 1, y = X2, col = name)) +
#   geom_boxplot() +
#   facet_grid(sim ~ n, scales = "free")

# Save plots -------------------------------------------------------------------
save(p_srs_a, p_srs_b, p_srs_c,
     p_srs_d, p_srs_e, p_srs_f,
     p_complex_a, p_complex_b, p_complex_c,
     p_complex_d, p_complex_e, p_complex_f,
     srs_plot, complex_plot,
     p_hist_a, p_hist_b, p_hist_c, p_hist_d,
     file = "vignettes/articles/simplots.RData")
