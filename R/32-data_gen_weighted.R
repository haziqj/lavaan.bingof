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

  population <-
    gen_data_bin(model_no = model_no, n = n * 25, seed = 123, H1 = H1,
                 return_all = TRUE)

  set.seed(seed)

  population <-
    population %>%
    # mutate(prob = 2 * rowSums(across(starts_with("eta"))),
    #        prob = exp(prob) / (1 + exp(prob)))
    mutate(prob = 1 / (1 + exp(.data$ystar1))) %>%
    select(-dplyr::contains("ystar"))

  sampled <-
    population %>%
    slice_sample(n = n, weight_by = prob, replace = FALSE) %>%
    # NOTE TO SELF: weight_by is actually probability of selection (sums to 1)
    mutate(wt = 1 / prob,
           wt = dplyr::n() * wt / sum(wt))

  if (isTRUE(return_all)) {
    sampled
  } else {
    sampled %>%
      select(wt, starts_with("y"))
  }
}
