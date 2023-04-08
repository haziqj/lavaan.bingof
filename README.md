
<!-- README.md is generated from README.Rmd. Please edit that file -->

# lavaan.bingof

<!-- badges: start -->

[![R-CMD-check](https://github.com/haziqj/lavaan.bingof/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/haziqj/lavaan.bingof/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

This is the accompanying R package for the paper

> Goodness-of-fit tests for composite likelihood estimation under simple
> random and complex survey sampling

## Usage

``` r
library(lavaan.bingof)

# Simulate binary data
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

# Fit lavaan model using PML estimation
(mod <- txt_mod(model.no = 1))
#> [1] "eta1 =~ NA*y1 + y2 + y3 + y4 + y5"
fit <- lavaan::sem(mod, dat, std.lv = TRUE, estimator = "PML")

# Test statistics
all_tests(fit)
#> # A tibble: 7 × 7
#>       i     W    df name           pval converged Omega2_rank
#>   <dbl> <dbl> <dbl> <chr>         <dbl> <lgl>           <int>
#> 1     1  4.96  5    Wald          0.420 TRUE               15
#> 2     1  1.81  3.14 WaldV2,MM3    0.638 TRUE               15
#> 3     1  4.96  5    WaldV3        0.421 TRUE               15
#> 4     1  4.07  4.53 Pearson       0.474 TRUE               15
#> 5     1  3.30  3.72 PearsonV2,MM3 0.465 TRUE               15
#> 6     1  4.11  4.34 RSS,MM3       0.441 TRUE               15
#> 7     1  4.93  5.00 Multn,MM3     0.424 TRUE               15
```
