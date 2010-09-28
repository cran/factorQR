\name{makeData}

\alias{makeData}

\title{Function to make synthetic data for the factorQR function}

\description{
  \code{makeData} simulates data from a factor quantile regression model.
}

\arguments{
  \item{N}{The sample size.}
  \item{whichFactor}{A vector that indicates which factor each manifest variable relates to.  E.g., \code{whichFactor = c(1,1,1,2,2)} would indicate a two-factor model, with the first three manifest variables relating to the first factor and the second two to the second factor.}
  \item{pQuant}{The quantile of interest. Defaults to \code{0.5}.}
  \item{lambda}{The vector of the non-zero elements of the factor loading matrix, with length equal to that of \code{whichFactor}.  Do not include the factor loadings related to the response variable.  Defaults to \code{1}.}
  \item{LambdaS}{The vector of factor loadings related to the response.  Must have length equal to the number of distinct values in \code{whichFactor}.  Defaults to \code{0}.}
  \item{Phi}{Matrix of latent factor covariances.  Must be symmetric and positive-definite and have dimension equal to the number of latent factors.  Defaults to the identity matrix.}
  \item{lapScale}{Scale of the asymmetric Laplace error distribution.  Defaults to \code{1}.}
  \item{Psi}{Vector of error variances for the manifest explanatory variables.}
}

\value{
  Returns a matrix whose first column is the response \code{Y} and whose remaining columns are the explanatory manifest variables with the underlying factor grouping implied by \code{whichFactor}.
}

\author{
  Lane F. Burgette, Department of Statistical Science,
  Duke University.  \email{lb131@stat.duke.edu}
}

