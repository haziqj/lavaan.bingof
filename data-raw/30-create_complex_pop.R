library(tidyverse)
library(lavaan.bingof)

the_seed <- 21324

# Under H0
SCHPOP_3f <- make_population(5, seed = the_seed, H1 = FALSE, Sigma2_attr = TRUE)
SCHPOP_1f8v <- make_population(2, seed = the_seed, H1 = FALSE, Sigma2_attr = TRUE)
SCHPOP_1f15v <- make_population(3, seed = the_seed, H1 = FALSE, Sigma2_attr = TRUE)

# Under H1 alternative hypothesis
SCHPOP_3f_alt <- make_population(5, seed = the_seed, H1 = TRUE, Sigma2_attr = TRUE)
SCHPOP_1f8v_alt <- make_population(2, seed = the_seed, H1 = TRUE, Sigma2_attr = TRUE)
SCHPOP_1f15v_alt <- make_population(3, seed = the_seed, H1 = TRUE, Sigma2_attr = TRUE)

# Save the populations
usethis::use_data(
  SCHPOP_3f, SCHPOP_1f8v, SCHPOP_1f15v,
  SCHPOP_3f_alt, SCHPOP_1f8v_alt, SCHPOP_1f15v_alt,
  overwrite = TRUE,
  compress = "xz"
)


