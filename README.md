
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lavaan.bingof

<!-- badges: start -->

[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![R-CMD-check](https://github.com/haziqj/lavaan.bingof/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/haziqj/lavaan.bingof/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/haziqj/lavaan.bingof/branch/main/graph/badge.svg)](https://app.codecov.io/gh/haziqj/lavaan.bingof?branch=main)
<!-- badges: end -->

<!-- https://github.com/r-lib/pkgdown/issues/133 -->

![](https://raw.githubusercontent.com/haziqj/lavaan.bingof/main/data-raw/mult_bern_data.png)

This is the accompanying R package for the research paper

> Jamil, H., Moustaki, I., & Skinner, C. (2023). Goodness-of-fit tests
> for composite likelihood estimation under simple random and complex
> survey sampling. *Manuscript in preparation*.

This package contains the functions to compute the test statistics and
conduct simulation studies described in the above manuscript. Currently,
the package implements the following tests based on univariate and
bivariate residuals of a binary factor analysis model:

|     | Name                    | R function          | Remarks                                          |
|-----|-------------------------|---------------------|--------------------------------------------------|
| 1   | Wald test               | `Wald_test()`       | Described in Reiser (1996)                       |
| 2   | Wald test (diagonal)    | `Wald_diag_test()`  | A more efficient Wald test                       |
| 3   | Wald test (VCOV free)   | `Wald_vcovf_test()` | Described in Maydeu-Olivares and Joe (2005,2006) |
| 4   | Pearson test            | `Pearson_RS_test()` | Rao-Scott adjusments                             |
| 5   | Pearson test            | `Pearson_test()`    | Moment matching approximation                    |
| 6   | Residual sum of squares | `RSS_test()`        | Moment matching approximation                    |
| 7   | Multinomial test        | `Multn_test()`      | Moment matching approximation                    |

## Installation

Install this package from this GitHub repository:

``` r
# install.packages("pak") 
pak::pkg_install("haziqj/lavaan.bingof")
library(lavaan.bingof)  # load package
```

Note: This package depends on the modified `{lavaan}` package version
0.6-14.9001, which can be installed from the GitHub repository at
[`haziqj/lavaan`](https://github.com/haziqj/lavaan).

## Usage

There are three main functionalities of this package:

1.  Generate simulated data either from an infinite population or from a
    finite population using a complex sampling procedure.

2.  Obtain the test statistic values, the degrees of freedom of these
    chi-square variates, and corresponding $p$-values to determine
    goodness-of-fit.

3.  Wrap functions 1 and 2 in a convenient way to perform simulation
    studies for Type I errors and power.

### Create a simulated data set of ordinal binary responses

``` r
(dat <- gen_data_bin(n = 1000, seed = 123))
#> # A tibble: 1,000 × 5
#>    y1    y2    y3    y4    y5   
#>    <ord> <ord> <ord> <ord> <ord>
#>  1 1     0     0     1     1    
#>  2 1     1     1     1     1    
#>  3 1     1     1     0     1    
#>  4 1     1     0     1     1    
#>  5 1     1     0     1     1    
#>  6 1     1     1     1     1    
#>  7 1     1     1     1     0    
#>  8 1     1     1     1     1    
#>  9 1     1     1     1     1    
#> 10 1     0     0     1     1    
#> # ℹ 990 more rows
```

### Obtain the various test statistics and $p$-values

``` r
# Fit lavaan model using PML estimation
(mod <- txt_mod(model_no = 1))
#> [1] "eta1 =~ NA*y1 + y2 + y3 + y4 + y5"
fit <- lavaan::sem(mod, dat, std.lv = TRUE, estimator = "PML")

# Test statistics
all_tests(fit)
#> # A tibble: 7 × 6
#>      X2    df name          pval W_rank     S
#>   <dbl> <dbl> <chr>        <dbl>  <int> <int>
#> 1 2.81   5    Wald         0.730     14    15
#> 2 0.862  3.31 WaldDiag,MM3 0.872     15    15
#> 3 2.80   5    WaldVCF      0.730      5    15
#> 4 2.46   4.41 PearsonRS    0.710     15    15
#> 5 1.86   3.63 Pearson,MM3  0.709     15    15
#> 6 2.30   4.18 RSS,MM3      0.707     15    15
#> 7 2.79   5.00 Multn,MM3    0.732     15    15
```

### Test statistics under a complex sampling scheme

``` r
# Simulate a two-stage stratified cluster sampling with 50 PSUs sampled per
# stratum, and 1 cluster sampled within each PSU.
(dat <- gen_data_bin_strcl(population = make_population(1), npsu = 50, 
                           seed = 9423))
#> # A tibble: 3,040 × 9
#>    type  school class    wt y1    y2    y3    y4    y5   
#>    <chr> <chr>  <chr> <dbl> <ord> <ord> <ord> <ord> <ord>
#>  1 A     A105   A105o 0.651 1     1     1     1     1    
#>  2 A     A105   A105o 0.651 1     1     1     1     1    
#>  3 A     A105   A105o 0.651 1     1     1     1     1    
#>  4 A     A105   A105o 0.651 1     1     1     1     1    
#>  5 A     A105   A105o 0.651 1     1     0     1     1    
#>  6 A     A105   A105o 0.651 1     1     1     1     1    
#>  7 A     A105   A105o 0.651 1     1     0     1     1    
#>  8 A     A105   A105o 0.651 1     1     1     1     1    
#>  9 A     A105   A105o 0.651 1     1     1     1     1    
#> 10 A     A105   A105o 0.651 1     1     1     1     1    
#> # ℹ 3,030 more rows

# Fit lavaan model and create survey object
fit0 <- lavaan::sem(mod, dat, std.lv = TRUE, estimator = "PML")  # ignore wt
fit1 <- lavaan::sem(mod, dat, std.lv = TRUE, estimator = "PML",
                    sampling.weights = "wt")
svy <- survey::svydesign(ids = ~ school + class, strata = ~ type,
                         weights = ~ wt, data = dat, nest = TRUE)

# Compare with and without sampling weights
Wald_test(fit0)
#>         X2 df name      pval W_rank  S
#> 1 4.561825  5 Wald 0.4716543     13 15
Wald_test(fit1, svy_design = svy)
#>         X2 df name      pval W_rank  S
#> 1 4.633965  5 Wald 0.4621618     13 15
```

### Simulation wrapper

``` r
# Conduct a simulation study based on a 5 factor model (32 repetitions only for
# illustration). Data generated according to a stratified complex sample.
(pc <- parallel::detectCores())   # how many cores do we have?
#> [1] 32

res <- run_ligof_sims(model_no = 1, nsim = pc, ncores = pc - 2, samp = "strat",
                      simtype = "type1")
#>|======================================================================| 100%
```

<!-- #>|======================================================================| 100% -->

``` r
res
#> 
#> ── LIGOF simulations ───────────────────────────────────────────────────────────
#> 
#> Settings
#> Number of replications: 32
#> Model: 1 (1F 5V)
#> Sampling design: Stratified sampling
#> Sample size: 1000
#> 
#> Simulations completed in 23.2 secs
```

``` r
summary(res)
#> 
#> ── LIGOF simulations summary ───────────────────────────────────────────────────
#> 
#> Model 1 (1F 5V) using stratified sampling design (n = 1000)
#> • Converged: 32 / 32
#> • Rank deficient: 0 / 32
#> • Significance level: 0.05
#> 
#> 
#> ============  ==============  =============  =========
#> Test name     Rejection rate  Mean X2 value  Mean d.f.
#> ============  ==============  =============  =========
#> Wald                   0.062           5.11       5.00
#> WaldDiag,MM3           0.000           2.99       3.44
#> WaldVCF                0.062           4.89       5.00
#> PearsonRS              0.062           3.75       4.15
#> Pearson,MM3            0.031           2.67       3.04
#> RSS,MM3                0.031           3.30       3.62
#> Multn,MM3              0.062           4.92       4.99
#> ============  ==============  =============  =========
```
