\encoding{latin1}
\name{care-package}
\alias{care-package}
\docType{package}
\title{The care Package}

\description{
The "care" package implements the CAR regression approach 
  described in Zuber and Strimmer (2011).  CAR scores measure the 
  correlations between the response and the Mahalanobis-decorrelated 
  predictors.  The squared CAR score is a natural measure of variable 
  importance and provides a canonical ordering of variables - see
  Zuber and Strimmer (2011) for details.
  
  This package provides functions for estimating CAR scores, for
  variable selection using CAR scores, and for estimating corresponding
  regression coefficients. Both shrinkage as well as empirical estimators 
  are available.

The name of the package refers to \bold{CAR} 
  \bold{e}stimation  and \bold{CAR} \bold{re}gression.
}

\author{Verena Zuber and Korbinian Strimmer (\url{https://strimmerlab.github.io/})}
\references{
  
Zuber, V., and K. Strimmer. 2011.  High-dimensional regression and 
variable selection using CAR scores.  Statist. Appl. Genet. Mol. Biol. 10: 34.
<DOI:10.2202/1544-6115.1730>

Website: \url{https://strimmerlab.github.io/software/care/}
}
\keyword{multivariate}
\seealso{
\code{\link{carscore}}, \code{\link{slm}}, 
 \code{\link{efron2004}}, \code{\link{lu2004}}.
}
