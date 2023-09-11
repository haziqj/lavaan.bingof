# Test statistics values (Stratified cluster sampling) Model 1

    Code
      res
    Output
      # A tibble: 6 x 6
           X2    df name          pval Xi_rank     S
        <dbl> <dbl> <chr>        <dbl>   <int> <int>
      1 1.82   5    Wald         0.874      14    15
      2 1.77   5    WaldVCF      0.880       5    15
      3 1.14   3.75 WaldDiag,MM3 0.865      15    15
      4 0.442  3.18 Pearson,MM3  0.944      15    15
      5 0.738  3.73 RSS,MM3      0.931      15    15
      6 1.75   4.99 Multn,MM3    0.881      15    15

# Test statistics values (Stratified cluster sampling) Model 2

    Code
      res
    Output
      # A tibble: 6 x 6
           X2    df name           pval Xi_rank     S
        <dbl> <dbl> <chr>         <dbl>   <int> <int>
      1  31.8 20    Wald         0.0454      34    36
      2  29.9 20    WaldVCF      0.0715      20    36
      3  12.5 10.5  WaldDiag,MM3 0.288       36    36
      4  11.5  9.51 Pearson,MM3  0.279       36    36
      5  16.4 12.3  RSS,MM3      0.190       36    36
      6  32.5 19.7  Multn,MM3    0.0352      36    36

# Test statistics values (Stratified cluster sampling) Model 4

    Code
      res
    Output
      # A tibble: 6 x 6
           X2    df name          pval Xi_rank     S
        <dbl> <dbl> <chr>        <dbl>   <int> <int>
      1 27.6   34   Wald         0.775      47    55
      2 26.4   34   WaldVCF      0.819      34    55
      3  7.28  13.8 WaldDiag,MM3 0.918      55    55
      4  7.17  11.4 Pearson,MM3  0.811      55    55
      5 12.1   18.0 RSS,MM3      0.841      55    55
      6 26.7   33.0 Multn,MM3    0.773      55    55

