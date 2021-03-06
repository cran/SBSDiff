% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sbs_diff.R
\name{sbs.log}
\alias{sbs.log}
\title{Satorra-Bentler Scaled Chi-Squared Difference Test (Based on Loglikelihood Values)}
\usage{
sbs.log(L0, L1, p0, p1, c0, c1)
}
\arguments{
\item{L0}{loglikelihood value for the more restrictive model (should be a negatige value)}

\item{L1}{loglikelihood value for the less restrictive model (should be a negatige value)}

\item{p0}{number of free parameters for the more restrictive model (with fewer freely estimated parameters)}

\item{p1}{number of free parametersfor the less restrictive model (with more freely estimated parameters)}

\item{c0}{scaling correction factor for the more restrictive model}

\item{c1}{scaling correction factor for the less restrictive model}
}
\value{
Change in model chi-squared, change in model degrees of freedom and the probability of rejecting the null
}
\description{
Takes loglikelihood values from nested models estimated using maximum likelihood
with robust standard errors, number of free parameters, scaling correlation factors and returns:
(1) Satorra-Bentler scaled change in model chi-squared (2) change in model degrees of freedom and
(3) the probability of rejecting the null.
}
\examples{
L0 <- -50

L1 <- -45

p0 <- 9

p1 <- 10

c0 <- 1

c1 <- 1

sbs.log(L0,L1,p0,p1,c0,c1)
}
