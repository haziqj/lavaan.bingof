# Test statistics values (Cluster sampling) Model 1

    Code
      res
    Output
      # A tibble: 6 x 6
           X2    df name          pval Xi_rank     S
        <dbl> <dbl> <chr>        <dbl>   <int> <int>
      1  3.53  5    Wald         0.618      14    15
      2  3.49  5    WaldVCF      0.626       5    15
      3  1.32  3.01 WaldDiag,MM3 0.725      15    15
      4  2.20  4.19 WaldDiag,RS2 0.726      15    15
      5  1.99  3.53 Pearson,MM3  0.667      15    15
      6  2.69  4.40 Pearson,RS2  0.671      15    15

# Test statistics values (Cluster sampling) Model 2

    Code
      res
    Output
      # A tibble: 6 x 6
           X2    df name          pval Xi_rank     S
        <dbl> <dbl> <chr>        <dbl>   <int> <int>
      1 15.2  20    Wald         0.763      34    36
      2 14.3  20    WaldVCF      0.812      20    36
      3  3.70  8.22 WaldDiag,MM3 0.895      36    36
      4  9.07 15.2  WaldDiag,RS2 0.883      36    36
      5  4.17 10.2  Pearson,MM3  0.947      36    36
      6  7.50 14.8  Pearson,RS2  0.937      36    36

# Test statistics values (Cluster sampling) Model 4

    Code
      res
    Output
      # A tibble: 6 x 6
           X2    df name          pval Xi_rank     S
        <dbl> <dbl> <chr>        <dbl>   <int> <int>
      1  42.5  34   Wald         0.150      50    55
      2  36.5  34   WaldVCF      0.354      34    55
      3  13.1  15.2 WaldDiag,MM3 0.612      55    55
      4  23.0  25.8 WaldDiag,RS2 0.621      55    55
      5  13.3  13.3 Pearson,MM3  0.448      55    55
      6  23.0  23.0 Pearson,RS2  0.460      55    55

