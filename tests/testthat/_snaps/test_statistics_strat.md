# Test statistics values (Stratified sampling) Model 1

    Code
      res
    Output
      # A tibble: 6 x 6
           X2    df name          pval Xi_rank     S
        <dbl> <dbl> <chr>        <dbl>   <int> <int>
      1 2.74   5    Wald         0.739      12    15
      2 2.74   5    WaldVCF      0.741       5    15
      3 0.664  3.04 WaldDiag,MM3 0.886      15    15
      4 1.52   3.90 Pearson,MM3  0.811      15    15
      5 1.94   4.33 RSS,MM3      0.789      15    15
      6 2.71   5.00 Multn,MM3    0.745      15    15

# Test statistics values (Stratified sampling) Model 2

    Code
      res
    Output
      # A tibble: 6 x 6
           X2    df name          pval Xi_rank     S
        <dbl> <dbl> <chr>        <dbl>   <int> <int>
      1 11.4  20    Wald         0.937      31    36
      2 11.2  20    WaldVCF      0.941      20    36
      3  3.04 10.1  WaldDiag,MM3 0.981      36    36
      4  9.06  8.53 Pearson,MM3  0.387      36    36
      5  9.59 11.9  RSS,MM3      0.644      36    36
      6 11.2  19.9  Multn,MM3    0.940      36    36

# Test statistics values (Stratified sampling) Model 4

    Code
      res
    Output
      # A tibble: 6 x 6
           X2    df name          pval Xi_rank     S
        <dbl> <dbl> <chr>        <dbl>   <int> <int>
      1  34.5 34    Wald         0.445      49    55
      2  31.9 34    WaldVCF      0.569      34    55
      3  14.0 13.1  WaldDiag,MM3 0.379      55    55
      4  12.2  9.98 Pearson,MM3  0.273      55    55
      5  18.3 15.7  RSS,MM3      0.285      55    55
      6  33.4 33.3  Multn,MM3    0.461      55    55

