% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_functions.R
\name{get_WAIC}
\alias{get_WAIC}
\title{Get Model Selection Statistics}
\usage{
get_WAIC(x, include_warning = FALSE)
}
\arguments{
\item{x}{\code{btpc_MCMC}, object output from performing MCMC using the \code{bTPC} function.}

\item{include_warning}{logical, should warnings be printed? Default is FALSE.}
}
\value{
A named vector containing the Watanabe-Akaike information criterion (WAIC), log point-wise predictive density (lppd), and pWAIC (a measure of the effective number of parameters being estimated) of the model.
}
\description{
Presents appropriate model selection in a convenient form.
}
\references{
Gelman, A., Hwang, J. & Vehtari, A. Understanding predictive information criteria for Bayesian models. Stat Comput 24, 997–1016 (2014). https://doi.org/10.1007/s11222-013-9416-2
}
