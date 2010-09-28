\name{checkFcn}

\alias{checkFcn}

\title{Check function}

\description{
  \code{checkFcn} is the check function, or tilted absolute
	value function.
}

\arguments{
  \item{x}{A vector of points at which we evaluate the function.}
  \item{p}{The quantile of interest.}
}

\value{
  Returns a vector with the same length as \code{x}.  The check loss is 
defined to be -1*(x < 0)*x*(1-p) + (x>0)*p*x.
}

\author{
  Lane F. Burgette, Department of Statistical Science,
  Duke University.  \email{lb131@stat.duke.edu}
}

