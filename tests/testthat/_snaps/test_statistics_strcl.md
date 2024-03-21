# Test statistics values (Stratified cluster sampling) Model 1

    Code
      res
    Output
      # A tibble: 6 x 6
           X2    df name          pval Xi_rank     S
        <dbl> <dbl> <chr>        <dbl>   <int> <int>
      1 1.82   5    Wald         0.873      14    15
      2 1.77   5    WaldVCF      0.880       5    15
      3 1.14   3.75 WaldDiag,MM3 0.864      15    15
      4 0.442  3.18 Pearson,MM3  0.944      15    15
      5 0.738  3.73 RSS,MM3      0.931      15    15
      6 1.76   4.99 Multn,MM3    0.881      15    15

# Test statistics values (Stratified cluster sampling) Model 2

    Code
      res
    Output
      # A tibble: 6 x 6
           X2    df name           pval Xi_rank     S
        <dbl> <dbl> <chr>         <dbl>   <int> <int>
      1  32.8 20    Wald         0.0355      34    36
      2  30.8 20    WaldVCF      0.0582      20    36
      3  12.6 10.5  WaldDiag,MM3 0.285       36    36
      4  11.6  9.51 Pearson,MM3  0.278       36    36
      5  16.4 12.3  RSS,MM3      0.189       36    36
      6  33.6 19.7  Multn,MM3    0.0266      36    36

# Test statistics values (Stratified cluster sampling) Model 4

    Code
      res
    Output
      # A tibble: 6 x 6
           X2    df name          pval Xi_rank     S
        <dbl> <dbl> <chr>        <dbl>   <int> <int>
      1 28.3   34   Wald         0.743      47    55
      2 27.1   34   WaldVCF      0.793      34    55
      3  7.30  13.8 WaldDiag,MM3 0.916      55    55
      4  7.18  11.4 Pearson,MM3  0.810      55    55
      5 12.1   18.0 RSS,MM3      0.841      55    55
      6 27.4   32.9 Multn,MM3    0.739      55    55

