dat <- gen_data_bin_wt()
fit <- lavaan::sem(model = txt_mod(1), data = dat, std.lv = TRUE, sampling.weights = "wt", estimator = "PML")

pp <- get_uni_bi_moments(fit)

pdot1 <-
  dat |>
  mutate(across(y1:y5, \(x) {
    wt * (x == 1)
  })) |>
  summarise(across(y1:y5, mean)) |>
  unlist()

pp$pdot1 - pdot1


Vy1 <- extract_lavaan_info(fit)$Var_ystar

Lambda <- matrix(partable(fit)$est[1:5], ncol = 1)
Theta  <- diag(partable(fit)$est[11:15])

Vy2 <- Lambda %*% t(Lambda) + Theta
Vy1



D0 <- get_Delta_mats(fit0)
D1 <- get_Delta_mats(fit1)


sim_X2 <- function(i) {
  dat <- gen_data_bin_wt()
  fit0 <- lavaan::sem(model = txt_mod(1), data = dat, std.lv = TRUE, estimator = "PML")
  fit1 <- lavaan::sem(model = txt_mod(1), data = dat, std.lv = TRUE, sampling.weights = "wt", estimator = "PML")

  extract_X2 <- function(fit, nparam = 10, wtd = TRUE) {

    tab <- lavaan::lavTables(fit, dimension = 0L) |>
      suppressWarnings()

    if (isTRUE(wtd)) {
      wtprop <-
        dat |>
        mutate(across(starts_with("y"), as.numeric)) |>
        unite("pattern", starts_with("y"), sep = "") |>
        summarise(wtprop = sum(wt), .by = pattern)
      sumw <- sum(wtprop$wtprop)
      wtprop$wtprop <- wtprop$wtprop / sumw

      tab <-
        left_join(tab, wtprop, by = "pattern") |>
          mutate(
            X2 = (wtprop - est.prop) ^ 2 / est.prop * sumw,
            G2 = 2 * wtprop * log(wtprop / est.prop) * sumw
          ) |>
        drop_na()

      # patternz <- expand.grid(rep(list(1:2), 5)) %>%
      #   apply(1, paste, collapse = "") %>%
      #   unique()
      #
      # identify_pattern <- function(x) {
      #   res <- rep(0, length(patternz))
      #   names(res) <- patternz
      #   res[patternz == x] <- 1
      #   res
      # }
      #
      # XX <-
      #   dat |>
      #   mutate(across(starts_with("y"), as.numeric)) |>
      #   unite("pattern", starts_with("y"), sep = "") |>
      #   mutate(id = row_number()) |>
      #   mutate(pp = map(pattern, identify_pattern)) |>
      #   unnest_wider(pp)
      #
      # Sigma <- cov.wt(select(XX, -wt, -pattern, -id), wt = XX$wt)$cov
      pihat <- tab$est.prop
      Sigma <- (diag(pihat) - tcrossprod(pihat)) %*% diag(1 / pihat)

      eval <- eigen(Sigma)$values
      dhat <- mean(eval)

      tab$X2 <- tab$X2 / dhat
      tab$G2 <- tab$G2 / dhat
    }

    c(
      X2 = sum(tab$X2),
      G2 = sum(tab$G2),
      df = 2 ^ 5 - 5 * (1 + 1) - 1
    )
  }

  tibble(
    wt = c("nowt", "wt"),
    X2 = list(
      extract_X2(fit0, wtd = FALSE),
      extract_X2(fit1, wtd = TRUE)
    )
  ) |>
    unnest_wider(X2) |>
    pivot_longer(c(X2, G2)) |>
    mutate(pval = pchisq(value, df, lower.tail = FALSE), i = i)
}

res <-
  future_map(
    1:1000,
    sim_X2,
    .options = furrr_options(seed = TRUE),
    .progress = TRUE
  )

res |>
  bind_rows() |>
  summarise(rej_rate = mean(pval < 0.05), .by = c(wt, name)) |>
  mutate(estimation = case_when(
    wt == "nowt" ~ "PML",
    wt == "wt" ~ "PMLW"
  )) |>
  select(-wt) |>
  pivot_wider(names_from = name, values_from = rej_rate)

# Check chi square test via LRT ------------------------------------------------
library(lavaan.bingof)
library(furrr)
plan(multisession, workers = parallel::detectCores() - 2)

sim_chisq <- function(i) {
  modno <- 1
  dat <- gen_data_bin_wt(model_no = modno, n = 5000)
  fit0 <- lavaan::sem(model = txt_mod(modno), data = dat, std.lv = TRUE, estimator = "PML")
  fit1 <- lavaan::sem(model = txt_mod(modno), data = dat, std.lv = TRUE, sampling.weights = "wt", estimator = "PML")


  make_chisq_tab <- function(fit) {
    tibble(
      thing = c("X2", "df"),
      X2 = as.numeric(fitMeasures(fit, c("chisq", "df"))),
      X2_scaled = as.numeric(fitMeasures(fit, c("chisq.scaled", "df.scaled")))
    ) |>
      pivot_longer(-thing) |>
      pivot_wider(id_cols = name, names_from = thing, values_from = value) |>
      mutate(pval = pchisq(X2, df, lower.tail = FALSE))
  }

  bind_cols(
    i = i,
    bind_rows(
      bind_cols(wt = "nowt", make_chisq_tab(fit0)),
      bind_cols(wt = "wt", make_chisq_tab(fit1))
    )
  )

}

res <-
  future_map(
    1:250,
    sim_chisq,
    .options = furrr_options(seed = TRUE),
    .progress = TRUE
  )

res |>
  bind_rows() |>
  summarise(
    rej_rate = mean(pval < 0.05),
    .by = c(wt, name)
  )
