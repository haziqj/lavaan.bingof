# Test statistics values (Cluster sampling) Model 1

    Code
      res
    Output
      # A tibble: 7 x 6
            X2    df name          pval Xi_rank     S
         <dbl> <dbl> <chr>        <dbl>   <int> <int>
      1 4.19    5    Wald         0.522      13    15
      2 1.17    2.97 WaldDiag,MM3 0.754      15    15
      3 2.60    5    WaldVCF      0.762       5    15
      4 0.762   3.27 PearsonRS    0.888      15    15
      5 0.0895  2.10 Pearson,MM3  0.963      15    15
      6 0.203   2.31 RSS,MM3      0.937      15    15
      7 1.66    3.53 Multn,MM3    0.734      15    15

# Test statistics values (Cluster sampling) Model 2

    Code
      res
    Output
      # A tibble: 7 x 6
            X2    df name             pval Xi_rank     S
         <dbl> <dbl> <chr>           <dbl>   <int> <int>
      1 164.   20    Wald         1.01e-24      33    36
      2   5.87  6.36 WaldDiag,MM3 4.81e- 1      36    36
      3  31.8  20    WaldVCF      4.51e- 2      20    36
      4   8.46  9.64 PearsonRS    5.50e- 1      36    36
      5   4.86  5.77 Pearson,MM3  5.33e- 1      36    36
      6   5.90  6.95 RSS,MM3      5.45e- 1      36    36
      7   3.08  2.45 Multn,MM3    2.86e- 1      36    36

# Test statistics values (Cluster sampling) Model 4

    Code
      res
    Output
      # A tibble: 7 x 6
            X2    df name             pval Xi_rank     S
         <dbl> <dbl> <chr>           <dbl>   <int> <int>
      1 491.   25    Wald         5.25e-88      45    55
      2   7.31  8.89 WaldDiag,MM3 5.94e- 1      55    55
      3 329.   34    WaldVCF      5.92e-50      34    55
      4  15.0  13.2  PearsonRS    3.19e- 1      55    55
      5   8.77  7.40 Pearson,MM3  3.06e- 1      55    55
      6   8.37  8.28 RSS,MM3      4.26e- 1      55    55
      7   9.22  2.05 Multn,MM3    1.05e- 2      46    55

