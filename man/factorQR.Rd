\name{factorQR}

\alias{factorQR}

\title{A Bayesian factor model for quantile regression}

\description{
  \code{factorQR} is an MCMC sampler to fit a Bayesian factor model for
 quantile regression.
}

\arguments{
  \item{factorForm}{A formula of the form \code{factorForm = Y ~ X1 + X2}, 
where \code{Y} is the response and variables on the right-hand side
are manifest variables related to the latent factors on which we are 
regressing \code{Y}.  The right-hand side variables will be centered automatically, though 
they are not scaled.}
  \item{nonFactorForm}{An optional formula of the form \code{nonFactorForm = ~ X3 + X4}. These covariates are used to model \code{Y}, but they do not relate to any of the latent factors.}
  \item{dataSet}{An optional data frame, list, or environment containing the variables in the model.}
  \item{pQuant}{Response quantile to model. Defaults to \code{pQuant=0.5}.}
  \item{whichFactor}{Vector of indicators to show factor grouping.  E.g., if \code{whichFactor=c(1,1,1,2,2,2)}, it would mean the first 
three variables in the right-hand side of \code{factorForm} group and the last three group together in two factors.  If \code{whichFactor = NULL}, a single latent factor will be assumed.}
  \item{nSamp}{Number of MCMC iterations, with a default of \code{5000}.}
  \item{burn}{Iterations of burn-in, with a default of \code{0}.}
  \item{thin}{Number of iterations to skip between stored values in the chain, with a default of \code{0}.}
  \item{cTau0}{Prior shape for \eqn{\tau}, which is the inverse scale of the response.  Defaults to \code{1}.}
  \item{dTau0}{Prior scale for \eqn{\tau}.}
  \item{cPsi0}{Prior shape for \eqn{\Psi}, the inverse scale of the manifest variables related to the factors.  Defaults to \code{1}.}
  \item{dPsi0}{Prior rate for \eqn{\Psi}. Defaults to \code{1}.}
  \item{sig0}{Hyperparameter for scale of free \eqn{\Lambda_{-s}} 
variables, which is rows of the factor loading matrix that do not 
correspond to the response.  Defaults to \code{1}.}
  \item{mu0}{Prior mean for \eqn{\Lambda_{-s}} components.  Defaults 
to \code{1}.}
  \item{R0}{Prior scale of \eqn{\Phi^{-1}}.  If not specified, defaults
to the identity matrix.}
  \item{nu0}{Prior degrees of freedom for inverse Wishart associated with \eqn{\Phi}.  Must be an integer greater than or equal to the dimension of \eqn{\Phi}, but defaults to one greater than that dimension.}
  \item{B0s}{Prior precision (i.e., inverse variance) for \eqn{\Lambda_{s}}.  Default is a diagonal matrix with non-zero values of \code{0.01}.  May be left at NULL, or changed to a non-negative scalar, a vector with 
length equal to the number of latent factors, or a symmetric, positive definite matrix with dimension equal to the number of latent factors.}
  \item{B0Beta}{Prior precision for regression coefficients related to \code{nonFactorForm}. Default is a diagonal matrix with non-zero values of \code{0.01}.  May be left at NULL, or changed to a non-negative scalar, a vector with length equal to the number of covariates not related to the latent factors, or a symmetric, positive definite matrix with that dimension. }
  \item{betaZero}{Starting value for \eqn{\beta}.}
  \item{PhiZero}{Initial value of \eqn{\Phi}.}
  \item{invPsiZero}{initial value of \eqn{\Psi^{-1}}.}
  \item{LambdaZero}{Initial value of \eqn{\Lambda_{-s}}. Must be a scalar or vector with length equal to the number of manifest variables in the right-hand side of \code{factorForm}.  Will be expanded to its matrix form internally, and elements that should be fixed at \code{1} will be adjusted if necessary, with warning.}
  \item{LambdaSZero}{Initial value of last row of \eqn{\Lambda}, i.e. the part that models the left-hand side of \code{factorForm}.  If specified, must either be a scalar or vector with length equal to the number of latent factors being modeled. Defaults to \code{0.1}}
  \item{OmegaZero}{Initial value for latent \eqn{\Omega}.}
  \item{verbose}{If \code{TRUE}, prints progress updates in Gibbs sampler.}
  \item{storeOmega}{If \code{TRUE}, stores the sampled \eqn{\Omega} values.}
  \item{latentInteract}{If \code{TRUE}, interacts the first two latent factors.  More latent/latent interactions are not currently supported.}
  \item{interactX}{If supplied, a matrix of variables with which to interact the latent factors.  Repeated columns are allowed.}
  \item{whichFactorInteract}{A vector that indicates which factors to interact with the columns of \code{interactX}.  For example, 
\code{whichFactorInteract = c(1,1)} would indicate that the two columns of \code{interactX} are to be interacted with the first latent factor.}
}

\value{
Returns an item of the class \code{factorQR} composed of the following components:
  \item{param}{Matrix of sampled parameter values.}
  \item{call}{The matched call.}
  \item{nReg}{The number of regression parameters.}
  \item{betLen}{The number of \eqn{\beta} components.}
  \item{nObs}{The number of observations.}
  \item{burn}{The number of Gibbs iterations before samples were stored.}
  \item{thin}{The number of Gibbs iterations between stored values.}
  \item{nSamp}{The total number of Gibbs iterations.}
  \item{nFact}{The number of modeled latent factors.}
  \item{nFactorX}{The number of manifest variables related to the factors.}
  \item{omega}{Sampled \eqn{\Omega} values, if \code{storeOmega} is \code{TRUE}.}
  \item{nFactInt}{Number of factor/manifest variable interactions.}
}

\author{
  Lane F. Burgette, Department of Statistical Science,
  Duke University.  \email{lb131@stat.duke.edu}
}

