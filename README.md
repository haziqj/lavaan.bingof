
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lavaan.bingof

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/haziqj/lavaan.bingof/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/haziqj/lavaan.bingof/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

![](analysis/mult_bern_data.png)

This is the accompanying R package for the paper

> Goodness-of-fit tests for composite likelihood estimation under simple
> random and complex survey sampling

This package contains the functions to compute the test statistics and
conduct simulation studies described in the above manuscript.

## Installation

Install this package from this GitHub repository:

``` r
# install.packages("pak") 
pak::pkg_install("haziqj/lavaan.bingof")
library(lavaan.bingof)  # load package
```

Note: This package depends on the modified `{lavaan}` package version
0.6-14.9001, which can be installed from the GitHub repository at
“[`haziqj/lavaan`](https://github.com/haziqj/lavaan)”.

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
#>  1 1     0     0     1     0    
#>  2 1     1     1     1     1    
#>  3 1     1     0     1     1    
#>  4 1     1     0     1     1    
#>  5 1     1     1     1     0    
#>  6 1     1     1     1     1    
#>  7 1     1     0     1     1    
#>  8 0     0     0     1     1    
#>  9 1     1     0     1     1    
#> 10 1     0     0     1     1    
#> # ℹ 990 more rows
```

### Obtain the various test statistics and $p$-values

``` r
# Fit lavaan model using PML estimation
(mod <- txt_mod(model.no = 1))
#> [1] "eta1 =~ NA*y1 + y2 + y3 + y4 + y5"
fit <- lavaan::sem(mod, dat, std.lv = TRUE, estimator = "PML")

# Test statistics
all_tests(fit)
#> # A tibble: 7 × 4
#>       W    df name           pval
#>   <dbl> <dbl> <chr>         <dbl>
#> 1  4.96  5    Wald          0.420
#> 2  1.81  3.14 WaldV2,MM3    0.638
#> 3  4.96  5    WaldV3        0.421
#> 4  4.07  4.53 Pearson       0.474
#> 5  3.30  3.72 PearsonV2,MM3 0.465
#> 6  4.11  4.34 RSS,MM3       0.441
#> 7  4.93  5.00 Multn,MM3     0.424
```

### Test statistics under a complex sampling scheme

``` r
# Simulate a two-stage stratified cluster sampling with 50 PSUs sampled per
# stratum, and 1 cluster sampled within each PSU.
(dat <- gen_data_bin_complex3(population = make_population(1), npsu = 50, 
                              seed = 9423))
#> # A tibble: 3,046 × 9
#>    type  school class    wt y1    y2    y3    y4    y5   
#>    <chr>  <int> <chr> <dbl> <ord> <ord> <ord> <ord> <ord>
#>  1 A          4 ae    0.796 1     1     1     1     1    
#>  2 A          4 ae    0.796 1     1     1     1     1    
#>  3 A          4 ae    0.796 1     1     0     0     1    
#>  4 A          4 ae    0.796 1     1     1     0     0    
#>  5 A          4 ae    0.796 1     1     1     1     1    
#>  6 A          4 ae    0.796 1     1     1     0     1    
#>  7 A          4 ae    0.796 1     1     1     1     1    
#>  8 A          4 ae    0.796 1     1     1     1     1    
#>  9 A          4 ae    0.796 1     1     1     1     1    
#> 10 A          4 ae    0.796 1     1     0     0     1    
#> # ℹ 3,036 more rows

# Fit lavaan model and create survey object
fit0 <- lavaan::sem(mod, dat, std.lv = TRUE, estimator = "PML")  # ignore wt
fit1 <- lavaan::sem(mod, dat, std.lv = TRUE, estimator = "PML",
                    sampling.weights = "wt")
svy <- survey::svydesign(ids = ~ school + class, strata = ~ type,
                         weights = ~ wt, data = dat, nest = TRUE)

# Compare with and without sampling weights
Wald_test(fit0)
#>          W df name      pval
#> 1 5.444331  5 Wald 0.3640888
Wald_test(fit1, svy_design = svy)
#>          W df name     pval
#> 1 6.194313  5 Wald 0.287768
```

### Simulation wrapper

``` r
# To conduct a simulation study based on a 5 factor model (25 repetitions only
# for illustration). Data generated according to a stratified complex sample.
(pc <- parallel::detectCores())   # how many cores do we have?
#> [1] 32
res <- ligof_sims(model.no = 1, nsim = pc, samp = "strat", simtype = "type1",
                  no.cores = pc)
#>|======================================================================| 100%
```

``` r
res %>%
  bind_rows() %>%
  group_by(name) %>%
  summarise(rej_rate = mean(pval < 0.1))
#> # A tibble: 7 × 2
#>   name          rej_rate
#>   <chr>            <dbl>
#> 1 Multn,MM3       0.0625
#> 2 Pearson         0.0312
#> 3 PearsonV2,MM3   0.0312
#> 4 RSS,MM3         0.0312
#> 5 Wald            0.0625
#> 6 WaldV2,MM3      0.0312
#> 7 WaldV3          0.0625
```
