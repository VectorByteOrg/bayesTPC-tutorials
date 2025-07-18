% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distribution_spec.R
\name{specify_likelihood}
\alias{specify_likelihood}
\title{Specify a Trait Data Likelihood}
\usage{
specify_likelihood(
  name = character(),
  formula = expression(),
  llh_parameters = character(),
  llh_constants = double()
)
}
\arguments{
\item{name}{character, the name of the new likelihood. Must not match any existing likelihood in bayesTPC.}

\item{formula}{expression, the likelihood code for trait data. The inferential parameter must be named 'm'.}

\item{llh_parameters}{optional named character, additional parameters used in the likelihood calculation (names) and their respective priors (value).}

\item{llh_constants}{optional named double, additional constants used in the likelihood calculation (names) and their respective default values (value).}
}
\value{
Returns a \code{btpc_likelihood} object, which can be modified or used in \code{\link[=specify_model]{specify_model()}} to create model specifications. Additionally, the name is registered so that the default specification can be referenced in the same fashion as the pre-included likelihoods.
}
\description{
Creates and registers a new likelihood for trait data. Useful for custom parameterizations of trait likelihoods, like using precision instead of variance in a normal distribution.
}
