# Figures for paper ------------------------------------------------------------

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
         sim = factor(sim, labels = c("1*F~5*V", "1*F~8*V", "1*F~15*V", "2*F~10*V", "3*F~15*V"))) %>%
  ggplot(aes(n, rej_rate, col = name, shape = name)) +
  geom_ribbon(aes(ymin = rej_rate - crit, ymax = rej_rate + crit,
                  fill = name), col = NA, alpha = 0.1) +
  geom_point() +
  geom_line() +
  geom_hline(aes(yintercept = alpha / 100), linetype = "dashed", col = "grey40") +
  scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 11, 16, 17, 15, 3, 7)) +
  facet_grid(sim ~ alpha_lab, labeller = label_parsed) +
  scale_x_continuous(breaks = unique(res_srs_power$n)) +
  scale_colour_viridis_d(option = "turbo", direction = -1) +
  scale_fill_viridis_d(option = "turbo", direction = -1) +
  guides(col = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
  labs(fill = NULL, col = NULL, shape = NULL,
       x = "Sample size (n)", y = "Rejection rate") +
  theme(legend.position = "bottom")

# ggsave("srs_plots.pdf", width = 7, height = 9)

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
         sim = factor(sim, labels = c("1*F~5*V", "1*F~8*V", "1*F~15*V", "2*F~10*V", "3*F~15*V"))) %>%
  ggplot(aes(n, rej_rate, col = name, shape = name)) +
  geom_ribbon(aes(ymin = rej_rate - crit, ymax = rej_rate + crit,
                  fill = name), col = NA, alpha = 0.1) +
  geom_point() +
  geom_line() +
  # geom_hline(aes(yintercept = alpha / 100), linetype = "dashed", col = "grey40") +
  scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 11, 16, 17, 15, 3, 7)) +
  facet_grid(sim ~ alpha_lab, labeller = label_parsed) +
  scale_x_continuous(breaks = unique(res_srs_power$n)) +
  scale_colour_viridis_d(option = "turbo", direction = -1) +
  scale_fill_viridis_d(option = "turbo", direction = -1) +
  guides(col = guide_legend(nrow = 1), shape = guide_legend(nrow = 1)) +
  labs(fill = NULL, col = NULL, shape = NULL,
       x = "Sample size (n)", y = "Rejection rate") +
  theme(legend.position = "bottom")

# ggsave("srs_plots_power.pdf", width = 7, height = 9)

# Complex plots
alpha <- 5
var_name <- paste0("rej_rate", alpha)
crit_name <- paste0("crit", alpha)

res_complex_type1 %>%
  filter(sampling != "Stratified") %>%
  mutate(sampling = factor(sampling, labels = c("Two-stage cluster", "Two-stage stratified cluster"))) %>%
  mutate(n = factor(n, labels = paste0("n =\n", unique(.data$n)))) %>%
  ggplot(aes(n, rej_rate5, fill = name,
             alpha = n_rank_def / .data$n_sims[1] * 100)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9) +
  geom_errorbar(aes(ymin = .data[[var_name]] - .data[[crit_name]],
                    ymax = .data[[var_name]] + .data[[crit_name]]),
                col = "black", width = 0.2, alpha = 1,
                position = position_dodge(width = 0.9)) +
  geom_hline(aes(yintercept = alpha / 100, linetype = "Nominal\nrej. level"),
             col = "grey50") +
  scale_linetype_manual(NULL, values = "dashed") +
  scale_alpha("% rank\ndef.", range = c(1, 0.5)) +
  theme(legend.position = "bottom") +
  labs(x = "Sample size", y = "Rejection proportion", fill = NULL,
       shape = NULL) +
  guides(fill = guide_legend(nrow = 2, order = 1),
         alpha = guide_legend(ncol = 3, order = 2)) +
  scale_fill_viridis_d(option = "turbo", direction = -1) +
  scale_colour_viridis_d(option = "turbo", direction = -1) +
  facet_grid(sim ~ sampling) +
  coord_cartesian(ylim = c(0, 0.1))
# ggsave("complex_type1.pdf", width = 7, height = 5)


res_complex_power %>%
  filter(sampling != "Stratified") %>%
  mutate(sampling = factor(sampling, labels = c("Two-stage cluster", "Two-stage stratified cluster"))) %>%
  mutate(n = factor(n, labels = paste0("n =\n", unique(.data$n)))) %>%
  ggplot(aes(n, rej_rate5, fill = name,
             alpha = n_rank_def / .data$n_sims[1] * 100)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.9) +
  geom_errorbar(aes(ymin = .data[[var_name]] - .data[[crit_name]],
                    ymax = .data[[var_name]] + .data[[crit_name]]),
                col = "black", width = 0.2, alpha = 1,
                position = position_dodge(width = 0.9)) +
  scale_linetype_manual(NULL, values = "dashed") +
  scale_alpha("% rank\ndef.", range = c(1, 0.5)) +
  theme(legend.position = "bottom") +
  labs(x = "Sample size", y = "Rejection proportion", fill = NULL,
       shape = NULL) +
  guides(fill = guide_legend(nrow = 2, order = 1),
         alpha = guide_legend(ncol = 3, order = 2)) +
  scale_fill_viridis_d(option = "turbo", direction = -1) +
  scale_colour_viridis_d(option = "turbo", direction = -1) +
  facet_grid(sim ~ sampling) +
  coord_cartesian(ylim = c(0, 1))
# ggsave("complex_power.pdf", width = 7, height = 5)


# Sankey plots -----------------------------------------------------------------
library(tidyverse)
library(ggsankey)
theme_set(theme_sankey())

dat <- gen_data_bin(model_no = 1, n = 1000, seed = 123) %>%
  mutate(across(starts_with("y"), function(x) as.numeric(x == 1)))

# Create bivariate moments
id <- combn(5, 2)
biv_dat <- NULL
for (k in seq_len(ncol(id))) {
  i <- id[1, k]
  j <- id[2, k]
  biv_dat[[k]] <- dat[[i]] * dat[[j]]
}
biv_dat <- as.data.frame(biv_dat)

# Create vector of name arrangements (for fill purposes as well)
uni_bi_y <- c(paste0("y", 5:1), rev(paste0("y", id[1, ], id[2, ])))
zy <-
  tibble(y = sort(uni_bi_y),
         z = rev(paste0("z", sprintf("%02d", 1:15)))) %>%
  arrange(z)
names(biv_dat) <-
  zy %>%
  slice(match(paste0("y", id[1, ], id[2, ]), y)) %>%
  pull(z)

get_uni2 <- function(x) {
  if (is.na(x)) return(NA)
  y <- zy %>%
    slice(match(x, z)) %>%
    pull(y)
  x <- paste0("y", str_split(gsub("y", "", y), "")[[1]])
  zy %>%
    slice(match(x, y)) %>%
    pull(z)
}

plot_df <-
  dat %>%
  mutate(pattern = apply(across(everything()), 1, paste0, collapse = ""),
         data = "All\n(n=1,000)") %>%
  bind_cols(biv_dat) %>%

  # Flow from response patterns to bivariate marginals (positive responses)
  pivot_longer(starts_with("z"), names_to = "bivariate",
               values_to = "bi_vals") %>%
  mutate(bivariate = case_when(bi_vals == 0 ~ NA, TRUE ~ bivariate)) %>%

  # Flow from bivariate to univariate positive marginals
  rowwise() %>%
  mutate(univariate = list(get_uni2(bivariate))) %>%
  unnest_longer(univariate) %>%

# Create the sankey data frame
  make_long(data, pattern, bivariate, univariate) %>%
  mutate(#col = my_col_fn(node),
    x = factor(x, labels = c("Multivariate\nBernoulli Data",
                             "Response\nPatterns", "Bivariate\nMoments",
                             "Univariate\nMoments")))

p_tmp <-
  ggplot(plot_df, aes(x = x, next_x = next_x, node = node,
                      next_node = next_node, fill = node)) +
  geom_sankey(flow.alpha = 0.5, space = 300, width = 0.25) +
  geom_sankey_text(aes(label = node), space = 300) +
  scale_fill_viridis_d(option = "turbo") +
  theme(legend.position = "none") +
  labs(x = NULL)

my_labeller <- function(x) {
  if (grepl("z", x)) {
    zy %>%
      slice(match(x, z)) %>%
      pull(y)
  } else {
    x
  }
}

# Extract data from plot builder
pb <- ggplot_build(p_tmp)
lab_dat <- pb$data[[3]] %>%
  as_tibble %>%
  rowwise() %>%
  mutate(labz = case_when(freq < 800 ~ NA, TRUE ~ my_labeller(node)))

p <-
  ggplot(plot_df, aes(x = x, next_x = next_x, node = node,
                      next_node = next_node,
                      fill = node)) +
  geom_sankey(flow.alpha = 0.5, space = 300, width = 0.3) +
  geom_label(data = lab_dat, aes(x, y, label = labz), fill = "gray80",
             alpha = 0.5, col = "gray20", inherit.aes = FALSE, size = 2.5) +
  scale_fill_viridis_d(option = "turbo") +
  scale_colour_viridis_d(option = "turbo") +
  theme(legend.position = "none") +
  labs(x = NULL, y = NULL); p

ggsave("mult_bern_data.png", p, width = 8.5, height = 5)
