#' Raw simulation results
#'
#' Raw simulation results
#'
#' For our simulation study, we varied the following settings:
#'
#' - Data generation mechanism: `srs`, `strat`, `clust`, `strcl`
#' - Sample size: 500, 1000, 2000, 3000
#' - Factor model: 5 models with various number of items and latent factors. See article details.
#' - Simulation type: Type I error analysis or Power analysis
#'
#' Thus, in total there are \eqn{4 \times 4 \times 5 \times 2 = 160} unique simulation scenarios. The R object `all_res` is a list of length 160 with each entry containing a [tibble()] of raw results (test statistic value, degrees of freedom, name of test, p-value, convergence, matrix ranks, and random seed) for each scenario. Note that there are 7 tests in total that was investigated, and each scenario was replicated a total of 1000 times.
#'
#' @name sim_results_all
#' @rdname sim_results_all
#' @format A list of length 160, each containing a large [tibble()] of raw results.
"all_res"

#' Summarised simulation results
#'
#' Summarised simulation results
#'
#' @name sim_results
#' @rdname sim_results
NULL

#' @rdname sim_results
#' @format A [tibble()] containing the summarised results of the simulation study segregated by test name (Wald, Wald V2, etc.), model (1F 5V, 1F 8V, etc.), and sample size.
#'
#'  - `n_sim` is the number of replications.
#'  - `n_converged` is the number of replications in which [lavaan::lavaan()] converged.
#'  - `n_rank_def` is the number of replications where the rank of \eqn{\Omega_2} matrix was found to be deficient.
#'  - `rej_rate_x` is the proportions of replications where the null hypothesis was rejected at the `x`% level.
#'
"res_srs_type1"

#' @rdname sim_results
#' @format NULL
"res_srs_power"

#' @rdname sim_results
#' @format NULL
"res_complex_type1"

#' @rdname sim_results
#' @format NULL
"res_complex_power"
