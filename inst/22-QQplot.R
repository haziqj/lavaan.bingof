library(tidyverse)
library(lavaan.bingof)
library(gt)

load("all_res.rda")

grab_sims <-
  function(
    x = all_res,
    samp = "srs",
    type = "type1",
    the_n = c(500, 1000, 2500, 5000, 10000)
  ) {

  mod_names <- c("1F 5V", "1F 8V", "1F 15V", "2F 10V", "3F 15V")
  type_of_sampling <- grepl(samp, names(x))
  type_of_analysis <- grepl(type, names(x))
  res <- list(NULL)

  for (n in the_n) {
    type_of_n <- grepl(paste0("n", n, "_"), names(x))
    ind <- which(type_of_sampling & type_of_analysis & type_of_n)
    tmp <- NULL
    for (i in seq_along(ind)) {
      tmp <- bind_rows(tmp, bind_cols(x[[ind[i]]], n = n, sim = mod_names[i]))
    }
    res <- c(res, list(tmp))
  }

  do.call("bind_rows", res) |>
    mutate(
      sim = factor(sim, levels = unique(sim)),
      name = factor(name, levels = unique(name)),
      alpha10 = pval < 0.1,
      alpha5 = pval < 0.05,
      alpha1 = pval < 0.01
    )
}

plot_X2_QQ <- function(samp = "srs") {
  dat <- grab_sims(samp = samp)
  ori_levs <- levels(dat$name)
  mod_levs <- gsub(",MM3", "", ori_levs)
  levels(dat$name) <- mod_levs

  dat |>
    mutate(
      df = mean(df),
      y = sort(X2),
      x = qchisq(ppoints(n()), df = first(df)),
      .by = c(name, sim, n)
    ) |>
    mutate(
      n = paste0("n = ", n) |>
        factor(levels = paste0("n = ", c(500, 1000, 2500, 5000, 10000)))
    ) |>
    ggplot(aes(x, y, col = n, shape = n)) +
    geom_point(alpha = 0.5) +
    geom_abline(col = "black") +
    facet_wrap(name ~ sim, scales = "free", ncol = 5) +
    theme_bw() +
    ggsci::scale_colour_d3() +
    labs(
      x = "Theoretical quantiles",
      y = "Observed quantiles",
      col = "Sample\nsize",
      shape = "Sample\nsize"
    ) +
    guides(
      override.aes = list(
        col = guide_legend(override.aes = list(alpha = 1))
      )
    )
}

create_tab_KS <- function(samp = "srs") {
  dat <- grab_sims(samp = samp)
  dat |>
    summarise(
      df = mean(df),
      KS = list(ks.test(X2, "pchisq", df = df,simulate.p.value = TRUE, B = 1000  )),
      # D = ks.test(X2, "pchisq", df = df)$statistic,
      # pval = ks.test(X2, "pchisq", df = df)$p.value,
      .by = c(name, sim, n)
    ) |>
    rowwise() |>
    mutate(
      D = KS$statistic,
      pval = KS$p.value,
      Dpval = glue::glue("{iprior::dec_plac(D, 3)} ({gtsummary::style_pvalue(pval)})")
    ) |>
    pivot_wider(id_cols = c(name, n),

      names_from = sim,
      values_from = Dpval
    ) |>
    arrange(name, n) |>
    gt(
      rowname_col = "n",
      groupname_col = "name"
    ) |>
    fmt_auto() |>
    tab_row_group(
      label = md("**Loadings**"),
      rows = contains("lambda")
    ) |>
    tab_options(quarto.use_bootstrap = FALSE)
}

# Save -------------------------------------------------------------------------

samps <- set_names(c("srs", "strat", "clust", "strcl"))
p_X2_QQ <- map(samps, plot_X2_QQ)
tab_KS <- map(samps, create_tab_KS)
save(p_X2_QQ, tab_KS, file = "qqplot_tab.RData")




