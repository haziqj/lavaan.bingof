library(tidyverse)
library(lavaan.bingof)


myfun <-     \(x, y, z, w) {

  cli::cli_alert_info("Now running {x} / {y} / pop_Sigma = {z} / Sigma2 = {w}")

  the_wt <- NULL
  if (y == "wt") the_wt <- "wt"

  out <- run_ligof_sims(
    model_no = 5,
    samp_size = 1000,
    nsim = 1000,
    samp = x,
    simtype = "type1",
    starting_seed = NULL,
    pop_Sigma = z,
    Sigma2 = w,
    the_wt = the_wt
  )
  cat("\n")
  summary(out)
}
myfun <- possibly(myfun, NA)

res <-
  expand_grid(
    samp = c("wtd", "strat", "clust", "strcl"),
    wt = c("wt", "nowt"),
    pop_Sigma = c(TRUE, FALSE),
    Sigma2 = c("theoretical", "weighted", "strat", "force_unweighted")
  ) |>
  filter(!(samp == "wtd" & pop_Sigma == TRUE)) |>
  filter(!(samp == "wtd" & Sigma2 == "strat")) |>
  filter(!(pop_Sigma == TRUE & Sigma2 != "theoretical")) |>
  filter(!(wt == "nowt" & Sigma2 != "theoretical")) |>
  filter(!(samp != "strat" & Sigma2 == "strat")) |>
  mutate(sim = pmap(
    list(samp, wt, pop_Sigma, Sigma2),
    myfun
  ))

save(res, file = "new_sims.RData")

# tmp <- myfun("strcl", "wt", FALSE, "weighted")
# res$sim[[20]] <- tmp

p_compare_sigma2 <-
  res |>
  mutate(
    sim = map(sim, ~ .x$tab)
  ) |>
  unnest(sim) |>
  mutate(
    Sigma2est = case_when(
      pop_Sigma == TRUE ~ "Sigma2 = Pop.",
      Sigma2 == "theoretical" ~ "Sigma2 = Theor.",
      Sigma2 == "weighted" ~ "Sigma2 = Wtd.",
      Sigma2 == "strat" ~ "Sigma2 = Strat.",
      Sigma2 == "force_unweighted" ~ "Sigma2 = Unwtd."
    ),
    Sigma2est = factor(Sigma2est, levels = c("Sigma2 = Theor.", "Sigma2 = Unwtd.", "Sigma2 = Wtd.", "Sigma2 = Strat.", "Sigma2 = Pop.")),
    wt = case_when(
      wt == "wt" ~ "Weighted",
      wt == "nowt" ~ "Ignore wt."
    ),
    samp = case_when(
      samp == "wtd" ~ "Informative",
      samp == "strat" ~ "Stratified",
      samp == "clust" ~ "Cluster",
      samp == "strcl" ~ "Strat. + Clust."
    ),
    samp = factor(samp, levels = c("Informative", "Stratified", "Cluster", "Strat. + Clust."))
  ) |>
  ggplot(aes(name, rej_rate, fill = name)) +
  geom_bar(stat = "identity", alpha = 0.8) +
  geom_hline(yintercept = 0.05, linetype = "dashed", col = "black", linewidth = 0.8) +
  facet_grid(samp * wt ~ Sigma2est) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  coord_cartesian(ylim = c(0, 0.2)) +
  labs(x = NULL, y = "Rejection rate", title = "Scenario 3F15V (n=1000): Rejection rates by simulation conditions", caption = "Total replications: 1000") +
  guides(fill = "none") +
  ggsci::scale_fill_jama()

save(p_compare_sigma2, file = "p_compare_sigma2.RData")
