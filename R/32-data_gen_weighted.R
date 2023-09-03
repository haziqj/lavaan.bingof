#' Generate an informative (weighted) data sample according to `model_no`
#'
#' @inheritParams gen_data_bin
#'
#' @return A [tibble()] containing ordinal binary values (0/1) for the items.
#' @export
#'
#' @examples
#' gen_data_bin_wt(1)
gen_data_bin_wt <- function(model_no = 1, seed = NULL, H1 = FALSE,
                            return_all = FALSE, n = 1000) {
  # if (!missing(n)) N <- n * 10
  res <-
    gen_data_bin(model_no = model_no, n = n * 10, seed = 1923, H1 = H1,
                 return_all = TRUE)

  set.seed(seed)

  res <-
    res %>%
    select(-dplyr::contains("ystar")) %>%
    mutate(prob = 2 * rowSums(across(starts_with("eta"))),
           prob = exp(prob) / (1 + exp(prob))) %>%
    slice_sample(n = n, weight_by = prob, replace = TRUE) %>%
    # NOTE TO SELF: weight_by is actually probability of selection (sums to 1)
    # mutate(sampled = stats::rbinom(dplyr::n(), size = 1, prob = prob) %>%
    #          as.logical()) %>%
    # filter(.data$sampled) %>%
    mutate(wt = 1 / prob,
           wt = dplyr::n() * wt / sum(wt))

  if (isTRUE(return_all)) {
    res
  } else {
    res %>%
      select(wt, starts_with("y"))
  }
}
