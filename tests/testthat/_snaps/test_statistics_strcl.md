# Test statistics values (Stratified cluster sampling) Model 1

    Code
      res
    Output
      # A tibble: 7 x 6
           X2    df name           pval W_rank     S
        <dbl> <dbl> <chr>         <dbl>  <int> <int>
      1 3.29   5    Wald          0.656     13    15
      2 0.871  2.66 WaldV2,MM3    0.782     15    15
      3 2.12   5    WaldV3        0.833      5    15
      4 0.957  3.51 Pearson       0.874     15    15
      5 0.298  2.41 PearsonV2,MM3 0.916     15    15
      6 0.457  2.69 RSS,MM3       0.899     15    15
      7 2.16   4.51 Multn,MM3     0.774     15    15

# Test statistics values (Stratified cluster sampling) Model 2

    Code
      res
    Output
      # A tibble: 7 x 6
            X2    df name              pval W_rank     S
         <dbl> <dbl> <chr>            <dbl>  <int> <int>
      1 214.   20    Wald          1.55e-34     31    36
      2   8.12  6.80 WaldV2,MM3    3.03e- 1     36    36
      3 110.   20    WaldV3        2.08e-14     20    36
      4  10.4   9.76 Pearson       3.85e- 1     36    36
      5   6.51  6.01 PearsonV2,MM3 3.69e- 1     36    36
      6   7.99  6.55 RSS,MM3       2.89e- 1     36    36
      7   3.11  1.86 Multn,MM3     1.90e- 1     36    36

# Test statistics values (Stratified cluster sampling) Model 4

    Code
      res
    Output
      # A tibble: 7 x 6
            X2    df name               pval W_rank     S
         <dbl> <dbl> <chr>             <dbl>  <int> <int>
      1 652.   27    Wald          3.43e-120     45    55
      2   5.19  7.85 WaldV2,MM3    7.23e-  1     55    55
      3 221.   34    WaldV3        2.60e- 29     34    55
      4   9.73 13.0  Pearson       7.15e-  1     55    55
      5   3.90  6.14 PearsonV2,MM3 7.06e-  1     55    55
      6   4.10  6.59 RSS,MM3       7.28e-  1     55    55
      7   1.69  1.36 Multn,MM3     2.78e-  1     48    55

