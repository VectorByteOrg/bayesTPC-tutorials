% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_spec.R
\name{specify_normal_model}
\alias{specify_normal_model}
\title{Specify model with normally distributed error}
\usage{
specify_normal_model(
  name = character(),
  parameters = character(),
  formula = expression(),
  constants = double(),
  dist_parameters = c(sigma.sq = "dexp(1)"),
  ...
)
}
\arguments{
\item{name}{character, The desired name of the model specification.}

\item{parameters}{named character, the names should correspond to the parameters being fit,
and the values should be the prior distributions to be drawn from for each respective parameter.
Uniform distributions should be used unless there is good reason to draw from another.}

\item{formula}{expression, The actual formula being fit.
Must include 'Temp' to represent temperature and all specified parameters and constants}

\item{constants}{optional named double, Represents any terms in the formula that should not be fit.}

\item{dist_parameters}{optional named character, additional parameters used in the likelihood calculation (names) and their respective priors (value).
Default value is 'c(sigma.sq = "dexp(1)")'.
If not provided, an exponential distribution with rate = 1 is used.}

\item{...}{Additional model specification attributes.}
}
\value{
Returns an object of type \code{btpc_normal_model}, which can then be used in other \code{bayesTPC} functions.
The model name is also registered, and so can be accessed using by passing only the name into functions.
However, user-defined models are not saved between sessions, and will be reset whenever the package is reloaded.
}
\description{
Creates an object with the required formatting to be fit using other \code{bayesTPC} functions.
}
\details{
\code{bayesTPC} does not verify if the priors specified are compatible with NIMBLE's dialect of BUGS.
All available distributions and formatting are provided on the
\href{https://r-nimble.org/html_manual/cha-writing-models.html#subsec:dists-and-functions}{NIMBLE user manual}.
}
\examples{
my_name <- "my_model"
my_formula <- expression(a * Temp^c + b)
my_parameters <- c(a = "dunif(0,1)", b = "dnorm(0,1)")
my_constants <- c(c = 1.5)
\dontrun{
my_model <- specify_normal_model(
  name = my_name,
  parameters = my_parameters,
  formula = my_formula,
  constants = my_constants
)
}
}
\seealso{
\code{\link[=remove_model]{remove_model()}}, \code{\link[=reset_models]{reset_models()}}, \code{\link[=specify_bernoulli_model]{specify_bernoulli_model()}}, \code{\link[=specify_model]{specify_model()}}, \code{\link[=specify_binomial_model]{specify_binomial_model()}}
}
