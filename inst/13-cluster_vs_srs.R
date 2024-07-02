library(tidyverse)
library(lavaan.bingof)
library(lme4)
library(furrr)
plan(multisession, workers = parallel::detectCores() - 1)

# Same population... 
# 1. SRS
# 2. cluster sample
# check power 

compare_srs_clust <- function(model_no = 1, samp_size = 1000, B = 1000) {
  
  cli::cli_alert_info("Running model {model_no} with sample size {samp_size}")
  pop <- lavaan.bingof:::make_population(model_no, H1 = TRUE)
  
  out <- furrr::future_map(
    .x = 1:B,
    .f = \(i) { 
      
      # Compare fit of 
      # fit0: SRS vs 
      # fit1: Cluster (unadjusted) vs 
      # fit2: Cluster (adjusted)
      
      dat_clust <- gen_data_bin_clust(pop, n = samp_size)
      fit1 <- lavaan::sem(
        model = txt_mod(model_no), 
        data = dat_clust, 
        estimator = "PML",
        std.lv = TRUE, 
        sampling.weights = NULL
      )
      fit2 <- lavaan::sem(
        model = txt_mod(model_no), 
        data = dat_clust, 
        estimator = "PML",
        std.lv = TRUE, 
        sampling.weights = "wt"
      )      
      fit0 <- lavaan::sem(
        model = txt_mod(model_no), 
        data = gen_data_bin(model_no, n = samp_size, H1 = TRUE),
        # data = gen_data_bin_srs(pop, n = samp_size),
        estimator = "PML",
        std.lv = TRUE
      )
      
      bind_rows(
        bind_cols(lavaan.bingof::all_tests(fit0, sim = i), samp = "srs"),
        bind_cols(lavaan.bingof::all_tests(fit1, sim = i), samp = "clust"),
        bind_cols(lavaan.bingof::all_tests(fit2, sim = i), samp = "clust-adj")
      ) 
    },
    .progress = TRUE,
    .options = furrr_options(seed = TRUE)
  )
  
  cat("\n")
  out
}
compare_srs_clust <- possibly(compare_srs_clust, NA)

res <-
  expand_grid(
    model_no = 1:5,
    samp_size = c(500, 1000, 2500, 5000, 10000)
  ) |>
  mutate(
    res = map2(model_no, samp_size, \(x, y) compare_srs_clust(x, y, B = 1000))
  )

save(res, file = "simulations/raw_power_cluster_vs_srs.RData")

res |>
  filter(sapply(res, \(x) !is.logical(x))) |>
  mutate(res = map(res, bind_rows)) |>
  unnest(res) |>
  summarise(
    rej_rate = mean(pval < 0.05, na.rm = TRUE),
    .by = c(model_no, samp_size, samp, name)
  ) |>
  pivot_wider(
    names_from = c(model_no, samp),
    values_from = rej_rate
  )

plot_df <- 
  res |> 
  filter(sapply(res, \(x) !is.logical(x))) |>
  mutate(res = map(res, bind_rows)) |>
  unnest(res) |>
  filter(converged) |>
  summarise(
    rej_rate = mean(pval < 0.05, na.rm = TRUE),
    .by = c(model_no, samp_size, samp, name)
  ) 
  # filter(name %in% c("Wald", "Pearson,MM3")) |>

tmp <- plot_df[plot_df$name == "WaldVCF" & plot_df$samp_size == 500 & plot_df$model_no %in% 3:5, ]
tmp$name <- "Wald"
plot_df[plot_df$name == "Wald" & plot_df$samp_size == 500 & plot_df$model_no %in% 3:5, ] <- tmp

tmp <- plot_df[plot_df$name == "WaldVCF" & plot_df$samp_size == 500 & plot_df$model_no %in% 3:5, ]
tmp$name <- "Wald"

plot_df$samp <- factor(plot_df$samp, levels = c("srs", "clust", "clust-adj"))
plot_df$model_no <- factor(plot_df$model_no, labels = c("1F5V", "1F8V",
                                                        "1F15V", "2F10V", "3F15V"))


p_srs_clust_compare <-
  plot_df |>
  filter(name %in% c("Wald", "Pearson,MM3")) |>
  ggplot(aes(samp_size, rej_rate, col = samp)) +
  geom_line(stat = "identity") +
  facet_grid(name ~ model_no) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(breaks = c(1000, 2500, 5000, 10000)) +
  labs(x = "Sample size", y = "Rejection rate", col = "Sampling\nmethod") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

save(p_srs_clust_compare, file = "simulations/p_srs_clust_compare.RData")
