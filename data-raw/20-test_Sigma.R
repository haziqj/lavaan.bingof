library(tidyverse)
library(survey)
data(api)
apiclus2 <- as_tibble(apiclus2)

# Two-stage cluster sampling specification:
#
# - Population is all schools in the state of California with > 100 students
# - Clusters are school districts (dname/dnum)
# - Elements are schools (sname/snum)
# - n = 40 districts SRS
# - One or more of schools in districts using SRS
# - fpc is the number of slements in population
api.twostage <- svydesign(id = ~dnum + snum, data = apiclus2,
                          fpc = ~fpc1 + fpc2)
summary(api.twostage)
