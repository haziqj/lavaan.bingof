# create_pairwise_table_tidyverse <- function(data, wt = NULL) {
#
#   p <- ncol(data)
#   if (is.null(wt)) data <- bind_cols(data, wt = 1)
#   else data <- bind_cols(data, wt = wt)
#
#   # Helper functions -----------------------------------------------------------
#   # pairwise_patterns <- function(x) {
#   #   expand_grid(yi = ordered(0:1),
#   #               yj = ordered(0:1)) %>%
#   #     unite(pattern, yi, yj, sep = "", remove = FALSE)
#   # }
#
#   count_patterns <- function(i, j, .data) {
#
#     selecting12wt <- function() dplyr::all_of(c(var1, var2, "wt"))
#     selecting1 <- function() dplyr::all_of(var1)
#     selecting2 <- function() dplyr::all_of(var2)
#
#     var1 <- paste0("y", i)
#     var2 <- paste0("y", j)
#     .data %>%
#       select(selecting12wt()) %>%
#       unite(pattern, -wt, sep = "", remove = FALSE) %>%
#       group_by(yi = dplyr::all_of(var1), yj = dplyr::all_of(var2), pattern) %>%
#       summarise(freq = n(), wtfreq = sum(wt), .groups = "drop") %>%
#       mutate(prop = freq / sum(freq), wtprop = wtfreq / sum(wtfreq))
#       # mutate(pattern = factor(pattern, levels = c("00", "01", "10", "11"))) %>%
#       # complete(yi, yj, pattern)
#   }
#
#   # The table ------------------------------------------------------------------
#   tab <-
#     as_tibble(t(combn(p, 2))) %>%
#     dplyr::rename("i" = 1, "j" = 2) %>%
#     mutate(combn = purrr::map2(i, j, \(x, y) count_patterns(x, y, .data = data) )) %>%
#     tidyr::unnest_wider(combn) %>%
#     tidyr::unnest_longer(c(yi, yj, pattern, freq, wtfreq, prop, wtprop))
#
#   S <- nrow(tab)
#
#   tab
# }
#
# create_pairwise_table_baseR <- function(data, wt = NULL) {
#
#   p <- ncol(data)
#
#   if (is.null(wt)) {
#     data$wt <- 1
#   } else {
#     data$wt <- wt
#   }
#
#   # Preallocate memory for the result
#   n_combinations <- choose(p, 2)
#   result_list <- vector("list", length = n_combinations)
#
#   k <- 1
#   for (i in 1:(p - 1)) {
#     for (j in (i + 1):p) {
#       var1 <- paste0("y", i)
#       var2 <- paste0("y", j)
#
#       # Using aggregate for weighted frequency counts
#       agg_df <- data.frame(yi = data[[var1]], yj = data[[var2]], freq = data$wt)
#       freq_df <- aggregate(freq ~ yi + yj, data = agg_df, length)
#       wtfreq_df <- aggregate(freq ~ yi + yj, data = agg_df, sum)
#       freq_df$wtfreq <- wtfreq_df$freq
#       freq_df$prop <- freq_df$freq / sum(freq_df$freq)
#       freq_df$wtprop <- freq_df$wtfreq / sum(freq_df$wtfreq)
#
#       # Store the result
#       result_list[[k]] <- list(i = i, j = j, table = freq_df)
#       k <- k + 1
#     }
#   }
#
#   # Convert result list to a data frame or other desired format
#   out <- lapply(result_list, as.data.frame) |> do.call(what = rbind)
#   names(out) <- gsub("table\\.", "", names(out))
#   out$pattern <- paste0(out$yi, out$yj)
#   as_tibble(out) %>%
#     select(i:yj, pattern, everything())
#
# }

create_pairwise_table_baseRv2 <- function(data, wt = NULL) {

  p <- ncol(data)

  if (is.null(wt)) {
    data$wt <- 1
  } else {
    data$wt <- wt
  }

  n_combinations <- choose(p, 2)
  result_list <- vector("list", length = n_combinations)

  # Pre-compute column names
  col_names <- paste0("y", 1:p)

  k <- 1
  for (i in 1:(p - 1)) {
    for (j in (i + 1):p) {
      var1 <- col_names[i]
      var2 <- col_names[j]

      # Create aggregate data frame
      agg_data <- data[c(var1, var2, "wt")]
      names(agg_data) <- c("yi", "yj", "freq")

      # Using tapply for aggregation
      freq <- tapply(agg_data$freq, list(agg_data$yi, agg_data$yj), length)
      wtfreq <- tapply(agg_data$freq, list(agg_data$yi, agg_data$yj), sum)

      # Convert to data frame and calculate proportions
      freq_df <- as.data.frame(as.table(freq))
      freq_df$wtfreq <- as.vector(wtfreq)
      freq_df$prop <- freq_df$Freq / sum(freq_df$Freq)
      freq_df$wtprop <- freq_df$wtfreq / sum(freq_df$wtfreq)

      # Store the result
      result_list[[k]] <- list(i = i, j = j, table = freq_df)
      k <- k + 1
    }
  }

  # Convert result list to a data frame or other desired format
  out <- do.call(rbind, lapply(result_list, as.data.frame))
  names(out) <- c("i", "j", "yi", "yj", "freq", "wtfreq", "prop", "wtprop")
  out$pattern <- paste0(out$yi, out$yj)
  out
}

# library(microbenchmark)
# microbenchmark(
#   create_pairwise_table_tidyverse(dat),
#   create_pairwise_table_baseR(dat),
#   create_pairwise_table_baseRv2(dat),
#   times = 10
# )

create_pairwise_table <- create_pairwise_table_baseRv2

pl_fn <- function(theta, model_no, data, wt = NULL) {

  # 1. Create pairwise table x
  # 2. Get Var(ystar) x
  # 3. Mutate pattern probabilities
  # 4. Output the pairwise likelihood

  convert_theta <- function(.theta) {
    lambdas <- .theta[grepl("lambda", names(.theta))]
    rhos <- .theta[grepl("rho", names(.theta))]
    taus <- .theta[grepl("tau", names(.theta))]

    list(lambdas = lambdas, rhos = rhos, taus = taus)
  }

  lambdas <- taus <- NULL  # to be overwritten by next line
  list2env(convert_theta(theta), environment())

  # Assign the parameters ------------------------------------------------------
  # loadings
  Lambda <- get_Lambda(model_no)
  Lambda[Lambda != 0] <- lambdas

  # factor correlations (if any)
  rhos <- 2 / (1 + exp(-rhos)) - 1  # expit [-1, 1]
  Psi <- cov_lv_mat(model_no)
  Psi[Psi != 1] <- rhos

  # thresholds
  tau <- taus

  # Var(ystar)
  neta <- ncol(Lambda)
  nitems <- nrow(Lambda)
  Theta <- matrix(0, nrow = nitems, ncol = nitems)
  diag(Theta) <- 1 - diag(Lambda %*% Psi %*% t(Lambda))
  Vy <- Lambda %*% Psi %*% t(Lambda) + Theta

  # Model probabilities --------------------------------------------------------
  calc_model_pairwise_prob <- function(i, j, pattern) {
    Vy_small <- Vy[c(i, j), c(i, j)]

    if (pattern == "00") {
      lower <- c(-Inf, -Inf)
      upper <- c(tau[i], tau[j])
    }
    if (pattern == "01") {
      lower <- c(-Inf, tau[j])
      upper <- c(tau[i], Inf)
    }
    if (pattern == "10") {
      lower <- c(tau[i], -Inf)
      upper <- c(Inf, tau[j])
    }
    if (pattern == "11") {
      lower <- c(tau[i], tau[j])
      upper <- c(Inf, Inf)
    }

    mnormt::sadmvn(lower = lower, upper = upper, mean = rep(0, 2),
                   varcov = Vy_small)
  }
  # res <-
  #   create_pairwise_table(data, wt) %>%
  #   mutate(prob = purrr::pmap_dbl(list(i, j, pattern), calc_model_pairwise_prob))
  res <- create_pairwise_table(data, wt)
  res$prob <- with(res, mapply(calc_model_pairwise_prob, i, j, pattern))


  # Pairwise log-likelihood ----------------------------------------------------
  # sum(res$wtfreq * log(res$prob))
  with(res, sum(wtfreq * log(wtprop / prob)))

}

get_Hinv_mat <- function(.data, .model_no, .wt = NULL, .fit = NULL) {
  theta0 <- get_true_values(.model_no)
  if (!is.null(.fit)) theta0[] <- coef(.fit)
  res <- optim(theta0, pl_fn, method = "BFGS", control = list(trace = 10),
               hessian = TRUE, model_no = .model_no, data = .data, wt = .wt)
  out <- nrow(.data)  * solve(res$hessian)
  attr(out, "coef") <- res$par
  out
}

# Check the functions ----------------------------------------------------------
# model_no <- 3
# samp_size <- 10000
# mod <- fit_facmod_pml(model_no, "strat", n = samp_size)
# pl_fn(get_true_values(model_no), model_no,
#       select(mod$dat, starts_with("y")), mod$dat$wt)
#
# library(profvis)
# profvis(pl_fn(get_true_values(model_no), model_no,
#               select(mod$dat, starts_with("y")), mod$dat$wt))
# profvis(create_pairwise_table(select(mod$dat, starts_with("y")), mod$dat$wt))

# Hinv_lav <- get_sensitivity_inv_mat(mod$fit)
# Hinv_pml <- get_Hinv_mat(select(mod$dat, starts_with("y")), model_no, mod$dat$wt)

# Does it improve the weighted test statistics? --------------------------------
#
#
#
# make_test_stats <- function(model_no = 1, n = 1000, pop) {
#   if(missing(pop)) pop <- make_population(model_no, seed = NULL)
#   # dat <- gen_data_bin_strat(pop, n = n, seed = NULL)
#   dat <- gen_data_bin_wt(model_no = model_no, n = n, seed = NULL)
#   fit <- sem(model = txt_mod(model_no), data = dat, estimator = "PML",
#              std.lv = TRUE, sampling.weights = "wt") %>%
#     suppressWarnings()
#   svy <- survey::svydesign(ids = ~ 1, weights = ~ wt, data = dat)
#
#   Hinv_pml <- lavaan.bingof:::get_Hinv_mat(select(dat, starts_with("y")),
#                                            model_no, dat$wt, fit)
#
#   all_tests(fit, svy, Hinv = Hinv_pml)
# }
# make_test_stats()
#
# res <- list()
# for (model_no in 1:2) {
#   cat(paste("\nRunning model", model_no, "\n"))
#   bigpop <- make_population(model_no)
#   plan(multisession, workers = 30)
#
#   possfn <- possibly(make_test_stats, NA)
#
#   out <-
#     future_map(1:250, ~possfn(model_no, n = 1000, pop = bigpop),
#                .progress = TRUE, .options = furrr_options(seed = NULL))
#
#   out <-
#     out[!sapply(out, is.na)] %>%
#     do.call(rbind, .) %>%
#     mutate(model_no = model_no)
#
#   res <- c(res, list(out))
# }
# res_df <- do.call(rbind, res)
# res_df %>%
#   group_by(model_no, name) %>%
#   summarise(rej_rate = mean(pval < 0.05)) %>%
#   pivot_wider(id_cols = name, names_from = model_no, values_from = rej_rate)








