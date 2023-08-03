# Test statistics values (Stratified cluster sampling) Model 1

    Code
      res
    Output
      # A tibble: 6 x 6
           X2    df name          pval Xi_rank     S
        <dbl> <dbl> <chr>        <dbl>   <int> <int>
      1  3.78  5    Wald         0.581      14    15
      2  3.72  5    WaldVCF      0.591       5    15
      3  2.80  3.37 WaldDiag,MM3 0.488      15    15
      4  3.83  4.48 WaldDiag,RS2 0.501      15    15
      5  2.57  3.18 Pearson,MM3  0.496      15    15
      6  3.53  4.24 Pearson,RS2  0.510      15    15

# Test statistics values (Stratified cluster sampling) Model 2

    Code
      res
    Output
      # A tibble: 6 x 6
           X2    df name           pval Xi_rank     S
        <dbl> <dbl> <chr>         <dbl>   <int> <int>
      1 29.8  20    Wald         0.0725      32    36
      2 25.5  20    WaldVCF      0.183       20    36
      3 17.1   9.58 WaldDiag,MM3 0.0611      36    36
      4 25.9  16.1  WaldDiag,RS2 0.0584      36    36
      5  7.25  7.81 Pearson,MM3  0.490       36    36
      6 12.4  13.1  Pearson,RS2  0.506       36    36

# Test statistics values (Stratified cluster sampling) Model 4

    Code
      res
    Output
      # A tibble: 6 x 6
           X2    df name              pval Xi_rank     S
        <dbl> <dbl> <chr>            <dbl>   <int> <int>
      1  74.0  34   Wald         0.0000867      50    55
      2  67.5  34   WaldVCF      0.000541       34    55
      3  24.4  14.7 WaldDiag,MM3 0.0533         55    55
      4  37.7  25.0 WaldDiag,RS2 0.0506         55    55
      5  21.6  12.0 Pearson,MM3  0.0422         55    55
      6  34.5  21.6 Pearson,RS2  0.0388         55    55

