library(tidyverse)
theme_set(theme_bw())

# If not downloaded, go to https://osf.io/2d97y/
load("res_srs_type1.rda")
load("res_srs_power.rda")
load("res_complex_type1.rda")
load("res_complex_power.rda")
load("simplots.RData")

# SRS Type I results
bind_rows(
  select(res_srs_type1, name, sim, n, rej_rate = rej_rate10, crit = crit10) %>%
    mutate(alpha = 10),
  select(res_srs_type1, name, sim, n, rej_rate = rej_rate5, crit = crit5) %>%
    mutate(alpha = 5),
  select(res_srs_type1, name, sim, n, rej_rate = rej_rate1, crit = crit1) %>%
    mutate(alpha = 1)
) %>%
  mutate(alpha_lab = factor(alpha, levels = c(10, 5, 1),
                            labels = c("alpha==10*symbol('%')",
                                       "alpha==5*symbol('%')",
                                       "alpha==1*symbol('%')")),
         name = factor(name, levels = rev(levels(name))),
         name = factor(name, labels = gsub(",MM3", "", levels(name))),
         name = factor(name, labels = gsub("Multn", "Multinomial", levels(name))),
         ok = !(!is.na(rej_rate) & !(abs(rej_rate - alpha / 100) < rej_rate)),
         N = factor(n),
         sim = factor(sim, labels = c("1*F~5*V", "1*F~8*V", "1*F~15*V",
                                      "2*F~10*V", "3*F~15*V"))) %>%
  ggplot(aes(rej_rate, name, shape = N, col = N)) +
  geom_vline(aes(xintercept = alpha / 100), linetype = "dashed") +
  geom_pointrange(aes(xmin = rej_rate - crit,
                      xmax = rej_rate + crit),
                  position = position_dodge(width = 0.5)) +
  facet_grid(sim ~ alpha_lab, labeller = label_parsed, scales = "free") +
  ggsci::scale_colour_d3() +
  # jcolors::scale_colour_jcolors(palette = "pal8") +
  # scale_colour_grey() +
  scale_alpha_manual(values = c(0.4, 1)) +
  scale_x_continuous(labels = scales::percent) +
  scale_shape_manual(values = c(16, 17, 15, 4, 1)) +
  labs(x = "Rejection rate", y = NULL, alpha = "Within\n95% interval",
       shape = "Sample size", col = "Sample size") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("srs_type1.pdf", width = 7, height = 9)

bind_rows(
  select(res_srs_power, name, sim, n, rej_rate = rej_rate10, crit = crit10) %>%
    mutate(alpha = 10),
  select(res_srs_power, name, sim, n, rej_rate = rej_rate5, crit = crit5) %>%
    mutate(alpha = 5),
  select(res_srs_power, name, sim, n, rej_rate = rej_rate1, crit = crit1) %>%
    mutate(alpha = 1)
) %>%
  mutate(alpha_lab = factor(alpha, levels = c(10, 5, 1),
                            labels = c("alpha==10*symbol('%')",
                                       "alpha==5*symbol('%')",
                                       "alpha==1*symbol('%')")),
         name = factor(name, labels = gsub(",MM3", "", levels(name))),
         name = factor(name, labels = gsub("Multn", "Multinomial", levels(name))),
         sim = factor(sim, labels = c("1*F~5*V", "1*F~8*V", "1*F~15*V", "2*F~10*V", "3*F~15*V"))) %>%
  ggplot(aes(n, rej_rate, col = name, shape = name)) +
  # geom_ribbon(aes(ymin = rej_rate - crit, ymax = rej_rate + crit,
  #                 fill = name), col = NA, alpha = 0.1) +
  geom_point() +
  geom_line(linewidth = 0.7) +
  # geom_hline(aes(yintercept = alpha / 100), linetype = "dashed", col = "grey40") +
  scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 11, 16, 17, 15, 3, 7)) +
  facet_grid(sim ~ alpha_lab, labeller = label_parsed) +
  scale_x_continuous(labels = scales::comma) +
  scale_colour_viridis_d(option = "turbo", direction = -1) +
  # scale_colour_grey() +
  scale_fill_viridis_d(option = "turbo", direction = -1) +
  guides(col = guide_legend(nrow = 2), shape = guide_legend(nrow = 2)) +
  labs(fill = NULL, col = NULL, shape = NULL, linetype = NULL,
       x = "Sample size (n)", y = "Rejection rate") +
  theme_bw() +
  theme(legend.position = "bottom", #panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        legend.key.width = unit(1.5,"cm"))

ggsave("srs_power.pdf", width = 7, height = 9)

# Complex plots
alpha <- 5
var_name <- paste0("rej_rate", alpha)
crit_name <- paste0("crit", alpha)

res_complex_type1 %>%
  filter(sampling != "Stratified") %>%
  mutate(sampling = factor(sampling, labels = c("Two-stage cluster",
                                                "Two-stage stratified cluster"))) %>%
  mutate(n = factor(n, labels = paste0("n =\n", unique(.data$n)))) %>%
  ggplot(aes(n, rej_rate5, fill = name,
             # alpha = n_rank_def / .data$n_sims[1] * 100
             )) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9, alpha = 0.8) +
  geom_errorbar(aes(ymin = .data[[var_name]] - .data[[crit_name]],
                    ymax = .data[[var_name]] + .data[[crit_name]]),
                col = "black", width = 0.2, alpha = 1,
                position = position_dodge(width = 0.9)) +
  geom_hline(aes(yintercept = alpha / 100, linetype = "Nominal\nrej. level")) +
  scale_linetype_manual(NULL, values = "dashed") +
  scale_alpha("% rank\ndef.", range = c(1, 0.5)) +
  theme(legend.position = "bottom") +
  labs(x = NULL, y = "Rejection proportion", fill = NULL,
       shape = NULL, title = "Type I errors") +
  guides(fill = guide_legend(nrow = 2, order = 1),
         alpha = guide_legend(ncol = 3, order = 2)) +
  scale_fill_viridis_d(option = "turbo", direction = -1) +
  # scale_fill_grey() +
  facet_grid(sim ~ sampling) +
  coord_cartesian(ylim = c(0, 0.15)) +
  theme_bw() +
  theme(legend.position = "none") -> p1; p1
# ggsave("complex_type1.pdf", width = 7, height = 5)

res_complex_power %>%
  filter(sampling != "Stratified") %>%
  mutate(sampling = factor(sampling, labels = c("Two-stage cluster",
                                                "Two-stage stratified cluster")),
         name = factor(name, labels = gsub(",MM3", "", levels(name))),
         name = factor(name, labels = gsub("Multn", "Multinomial", levels(name)))) %>%
  mutate(n = factor(n, labels = paste0("n =\n", unique(.data$n)))) %>%
  ggplot(aes(n, rej_rate5, fill = name,
             #alpha = n_rank_def / .data$n_sims[1] * 100
             )) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9, alpha = 0.8) +
  geom_errorbar(aes(ymin = .data[[var_name]] - .data[[crit_name]],
                    ymax = .data[[var_name]] + .data[[crit_name]]),
                col = "black", width = 0.2, alpha = 1,
                position = position_dodge(width = 0.9)) +
  scale_linetype_manual(NULL, values = "dashed") +
  scale_alpha("% rank\ndef.", range = c(1, 0.5)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  labs(x = "Sample size", y = "Rejection proportion", fill = NULL,
       shape = NULL, title = "Power analysis") +
  guides(fill = guide_legend(nrow = 1, order = 1),
         alpha = guide_legend(ncol = 3, order = 2)) +
  scale_fill_viridis_d(option = "turbo", direction = -1) +
  # scale_fill_grey() +
  facet_grid(sim ~ sampling) +
  coord_cartesian(ylim = c(0, 1)) -> p2; p2

# Save plots -------------------------------------------------------------------
cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(0.8, 1))
ggsave("complex_plot.pdf", width = 7, height = 9)
ggsave("complex_type1.png", p1, width = 7, height = 5)
ggsave("complex_power.png", p2, width = 7, height = 5)
# ggsave("hist_1_srs.pdf", p_hist_a + theme(panel.grid = element_blank()),
#        width = 7, height = 9)
# ggsave("hist_1_strat.pdf", p_hist_b + theme(panel.grid = element_blank()),
#        width = 7, height = 9)
# ggsave("hist_1_clust.pdf", p_hist_c + theme(panel.grid = element_blank()),
#        width = 7, height = 9)
# ggsave("hist_1_strcl.pdf", p_hist_d + theme(panel.grid = element_blank()),
#        width = 7, height = 9)

