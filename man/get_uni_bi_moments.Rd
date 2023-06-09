% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/40-test_statistics.R
\name{get_uni_bi_moments}
\alias{get_uni_bi_moments}
\title{Get univariate and bivariate moments}
\usage{
get_uni_bi_moments(.lavobject)
}
\arguments{
\item{.lavobject}{A \code{\link[lavaan:lavaan]{lavaan::lavaan()}} fit object.}
}
\value{
A list of univariate and bivariate moments.
}
\description{
Returns univariate and bivariate moments (i.e. positive
probabilities only) based on model i.e. \code{pidot1} and \code{pidot2} and
(weighted) sample i.e. \code{pdot1} and \code{pdot2}.
}
\examples{
fit <- lavaan::sem(txt_mod(1), gen_data_bin(1, n = 500), std.lv = TRUE,
                   estimator = "PML")
get_uni_bi_moments(fit)
}
\seealso{
\code{\link[=get_theoretical_uni_bi_moments]{get_theoretical_uni_bi_moments()}}
}
