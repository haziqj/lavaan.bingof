# Test statistics values (Stratified sampling) Model 1

    Code
      res
    Output
      # A tibble: 6 x 6
           X2    df name          pval Xi_rank     S
        <dbl> <dbl> <chr>        <dbl>   <int> <int>
      1 2.75   5    Wald         0.738      12    15
      2 2.74   5    WaldVCF      0.739       5    15
      3 0.665  3.04 WaldDiag,MM3 0.886      15    15
      4 1.52   3.90 Pearson,MM3  0.811      15    15
      5 1.94   4.33 RSS,MM3      0.789      15    15
      6 2.71   5.00 Multn,MM3    0.744      15    15

# Test statistics values (Stratified sampling) Model 2

    Code
      res
    Output
      # A tibble: 6 x 6
           X2    df name          pval Xi_rank     S
        <dbl> <dbl> <chr>        <dbl>   <int> <int>
      1 11.5  20    Wald         0.933      31    36
      2 11.3  20    WaldVCF      0.938      20    36
      3  3.06 10.1  WaldDiag,MM3 0.981      36    36
      4  9.11  8.56 Pearson,MM3  0.386      36    36
      5  9.62 11.9  RSS,MM3      0.644      36    36
      6 11.3  19.9  Multn,MM3    0.936      36    36

# Test statistics values (Stratified sampling) Model 4

    Code
      res
    Output
      # A tibble: 6 x 6
           X2    df name          pval Xi_rank     S
        <dbl> <dbl> <chr>        <dbl>   <int> <int>
      1  35.7 34    Wald         0.388      49    55
      2  33.0 34    WaldVCF      0.517      34    55
      3  14.1 13.1  WaldDiag,MM3 0.376      55    55
      4  12.2  9.96 Pearson,MM3  0.272      55    55
      5  18.3 15.7  RSS,MM3      0.284      55    55
      6  34.7 33.3  Multn,MM3    0.401      55    55

