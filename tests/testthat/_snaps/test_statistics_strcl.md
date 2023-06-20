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
      4 1.69   4.58 WaldDiag,RS2 0.855      15    15
      5 0.442  3.18 Pearson,MM3  0.944      15    15
      6 1.04   4.19 Pearson,RS2  0.917      15    15

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
      4  18.8 16.3  WaldDiag,RS2 0.296       36    36
      5  11.5  9.51 Pearson,MM3  0.279       36    36
      6  16.7 14.2  Pearson,RS2  0.286       36    36

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
      4 15.8   24.5 WaldDiag,RS2 0.909      55    55
      5  7.17  11.4 Pearson,MM3  0.811      55    55
      6 15.6   21.4 Pearson,RS2  0.809      55    55

