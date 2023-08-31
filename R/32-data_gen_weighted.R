#' Generate an informative (weighted) data sample according to `model_no`
#'
#' @inheritParams gen_data_bin
#' @param N The population size to sample from.
#' @param n The sample size, which is a subsample of `N`. If this is supplied,
#'   then `N` will be ignored.
#'
#' @return A [tibble()] containing ordinal binary values (0/1) for the items.
#' @export
#'
#' @examples
#' gen_data_bin_wt(1)
gen_data_bin_wt <- function(model_no = 1, seed = 123, H1 = FALSE,
                            return_all = FALSE, N = 1000, n) {

  if (!missing(n)) N <- n * 2

  res <-
    gen_data_bin(model_no = model_no, n = N, seed = seed, H1 = H1,
                 return_all = TRUE) %>%
    select(-dplyr::contains("ystar")) %>%
    mutate(prob = rowSums(across(starts_with("eta"))),
           prob = 1 / (1 + exp(prob))) %>%
    mutate(sampled = stats::rbinom(dplyr::n(), size = 1, prob = prob) %>%
             as.logical()) %>%
    filter(.data$sampled) %>%
    mutate(wt = 1 / prob,
           wt = dplyr::n() * wt / sum(wt))

  if (isTRUE(return_all)) {
    res
  } else {
    res %>%
      select(wt, starts_with("y"))
  }
}

