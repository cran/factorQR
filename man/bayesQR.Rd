\name{bayesQR}

\alias{bayesQR}

\title{Bayesian quantile regression}

\description{
  \code{bayesQR} is an MCMC sampler to fit a Bayesian quantile regression model.  This does not assume a factor structure.
}

\arguments{
  \item{formula}{A formula of the form \code{formula = Y ~ X1 + X2}, 
where \code{Y} is the response and variables on the right-hand side
are covariates.}
  \item{dataSet}{An optional data frame, list, or environment containing the variables in the model.}
  \item{pQuant}{Response quantile to model. Defaults to \code{pQuant=0.5}.}
  \item{nSamp}{Number of MCMC iterations, with a default of \code{5000}.}
  \item{burn}{Iterations of burn-in, with a default of \code{0}.}
  \item{thin}{Number of iterations to skip between stored values, with a default of \code{0}.}
  \item{C0}{Prior shape for \eqn{\tau}, which is the inverse scale of the response.  Defaults to \code{1}.}
  \item{D0}{Prior scale for \eqn{\tau}.}
  \item{B0}{Prior precision (i.e., inverse variance) for \eqn{\beta} regression parameters.  Default is a diagonal matrix with non-zero values of \code{0.01}.  May be left at NULL, or changed to a non-negative scalar, a vector with 
length equal to the number of covariates, or a symmetric, positive semi-definite matrix with dimension equal to the number of covariates.}
  \item{betaZero}{Starting value for \eqn{\beta}.}
  \item{verbose}{If \code{TRUE}, prints progress updates in Gibbs sampler.}
}

\value{
Returns an item of the class \code{bayesQR} composed of the following components:
  \item{param}{Matrix of sampled parameter values.}
  \item{call}{The matched call.}
  \item{betLen}{The number of \eqn{\beta} components.}
  \item{nObs}{The number of observations.}
  \item{burn}{The number of Gibbs iterations before samples were stored.}
  \item{thin}{The number of Gibbs iterations between stored values.}
  \item{nSamp}{The total number of Gibbs iterations.}
}

\author{
  Lane F. Burgette, Department of Statistical Science,
  Duke University.  \email{lb131@stat.duke.edu}
}

