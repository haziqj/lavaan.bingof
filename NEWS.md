# lavaan.bingof (development version)

## 2023-09-10

- Just as a sanity check, I coded my own pairwise likelihood function to see if the weighted PML in `{lavaan}` is doing what is expected. Short answer is yes--no issues there, including the $H^{-1}$ matrix. However no improvement to test statistics.
- Another idea was to fix the parameter values in `{lavaan}`. Unfortunately the fact that the model has no free parameters means that no Hessian and derivatives are calculated.
- Tested whether the option `sampling.weights.normalization = "none"` had any effect on test statistics. The answer is no. Exact same values are obtained.

## 2023-06-04

- Rewrote the function for estimating complex multinomial matrix $\boldsymbol\Sigma$ which uses the method from the `{survey}` package. (Although the one from that package does not work well when $p=15$ and $S=120$).
- Added bootstrap functionality for finding the complex multinomial matrix. This uses the bootstrap sampling functionality from the `{svrep}` package.
- Increment package version to 0.1.1

## 2023-05-12

- Bug in the `gen_data_bin_strat()` function. It ignores the sample size argument `n` so all simulations were actually using the same sample size.
- Update simulation results accordingly.

## 2023-04-15

- Added test suite for the package.
- Change of notation: The test statistic (previously $W$) is now denoted $X^2$. This makes sense because they are all generally chi-square variates. Avoids confusion with 'weights'.

## 2023-04-14

- Fixed bug in `lavaan.bingof:::create_Sigma2_matrix_complex()` where it did not loop over strata.
- Fixed bug in `lavaan.bingof:::extract_lavaan_info()` would extract the $2^p$ proportion tables. No wonder the simulations ran slow... oops.
- Now using `MASS::ginv` for inverting $\boldsymbol \Xi$ weights. Not sure if this is a good idea? This really only affects the Wald test and Multinomial test (apparently anyway).
- Re-ran all simulations 
- Wald and Wald V3 performs very poorly when the rank of weight matrix $\boldsymbol\Xi$ is rank deficient. Noticed that this occurs often in complex sampling and $\boldsymbol\Sigma_2$ matrix is not full rank. 

## 2023-04-08

- Create R package `{lavaan.bingof}` for easier shipment of R codes to implement the goodness-of-fit tests.
- Migrate from `{bookdown}` format presentation to `{pkgdown}` articles. 
- Previous repo [haziqj/gof-comp-lik](https://github.com/haziqj/gof-comp-lik/) and accompanying website will likely be shutdown.
 
## 2023-04-01

- After further discussions, decide to stratify the population based on the latent variable $\eta$. See technical details for clarification.
- In turn, the only missing component was the proper estimation of the multinomial matrix $\boldsymbol\Sigma$ using sampling weights.
- Finally, the complex sampling simulation works!

## 2023-03-22

- Experimenting with a different complex sampling scheme, whereby the population is determined by multiple variable including stratum-level, cluster-level, and other individual characteristics (called shape, colour, and "mystery").
- Inspired by the blog post "[Estimating relative risk in a simulated complex survey by @ellis2013nz](https://www.r-bloggers.com/2018/08/estimating-relative-risk-in-a-simulated-complex-survey-by-ellis2013nz/)".
- However, generating data this way makes the population highly misspecified, and so none of the test statistics work.

## 2023-02-26

- Wrote the technical details section for the derivation of the test statistics using both maximum likelihood estimation and pairwise maximum likelihood estimation.
- Consequently all R code was rewritten for efficiency.
- Added the missing parts of the test statistics, including
	- Degree of freedom calculation via Rao-Scott type adjustment and moment matching.
	- New tests: Wald V3, RSS and Multinomial.

## 2023-02-20

- Decided on a very basic stratified sampling procedure, where the population consists of 50 strata with varying population size within each stratum. The number of PSU to be selected is 20 so total sample size is 1000.
- To create homogeneity within stratum, a stratum level random effect was generated.
- Although bias was reduced when sampling weights were used, the test statistics values still were suspect.

## 2023-01-25

- After a long hiatus, revisit the sampling weights problem.
- Successfully modified the relevant parts in `{lavaan}` so that sampling weights can be used. The objective function (to be minimised) now takes the normalised weights instead of frequency counts.
- See github.com/haziqj/lavaan and install this fork of `{lavaan}`
- Next problem is how best to generate a complex sample?

## 2021-09-25

- Simulation results reproduced for the SRS case.
- Unable to proceed with complex sampling case because `{lavaan}` does not accept `sampling.weights` for `estimator = "PML"`.

## 2021-08-29

- Inherit R scripts and rewrite for reproducibility.
- All R code tested.
