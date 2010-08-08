\encoding{latin1}
\name{care-package}
\alias{care-package}
\docType{package}
\title{The care Package}

\description{
The "care" package implements the CAR regression approach described 
in Zuber and Strimmer (2010).
CAR scores measure the correlations between the response and the 
 Mahalanobis-decorrelated predictors.  The squared CAR score is a
natural measure of variable importance and provides
a canonical ordering of variables - see
\url{http://arxiv.org/abs/1007.5516} for details.

This package provides functions for estimating CAR scores, for
variable selection using CAR scores, and for estimating corresponding
regression coefficients. Both shrinkage as well as empirical estimators 
are available.

The name of the package refers to \bold{CAR} 
  \bold{e}stimation  and \bold{CAR} \bold{re}gression.
}

\author{Verena Zuber and Korbinian Strimmer (\url{http://strimmerlab.org/})}
\references{
  See website: \url{http://strimmerlab.org/software/care/}
}
\keyword{multivariate}
\seealso{
\code{\link{care}}, \code{\link{carscore}}, \code{\link{fitlm}}, \code{\link{efron2004}}, \code{\link{lu2004}}.
}
