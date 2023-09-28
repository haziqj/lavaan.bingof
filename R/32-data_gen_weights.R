# population should have 50 strata of unequal sizes Nh. Total pop size is roughly N = 1,000,000.
# strat sample from each strata, select say 20 individuals. The prob of selection is 20 / Nh


#' Title
#'
#' @param model_no
#' @param seed
#' @param H1
#' @param return_all
#' @param N
#' @param n
#'
#' @return
#' @export
#'
#' @examples
gen_data_bin_wt <- function(model_no = 1, seed = 123, H1 = FALSE,
                            return_all = FALSE, N = 1000, n) {

  if (!missing(n)) N <- n * 2

  gen_data_bin(model_no = model_no, n = N, seed = seed, H1 = H1,
                 return_all = TRUE) %>%
    select(-contains("ystar")) %>%
    mutate(prob = rowSums(across(starts_with("eta"))),
           prob = 1 / (1 + exp(prob))) %>%
    mutate(sampled = rbinom(dplyr::n(), size = 1, prob = prob) %>%
             as.logical()) %>%
    filter(sampled) %>%
    mutate(wt = 1 / prob,
           wt = dplyr::n() * wt / sum(wt)) %>%
    select(wt, starts_with("y"))

}

# dat <- gen_data_bin_wt()
# svy <- svydesign(ids = ~ 1, weights = ~ wt, data = dat)
