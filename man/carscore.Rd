\name{carscore}
\alias{carscore}

\title{Estimate CAR Scores and Marginal Correlations}
\usage{
carscore(x, y, estimator=c("empirical", "shrinkage"), diagonal=FALSE, verbose=TRUE)
}
\arguments{
  \item{x}{matrix of predictors (columns correspond to variables).}
  
  \item{y}{univariate response variable.}

  \item{estimator}{either \code{"empirical"} (the default) or \code{"shrinkage"}.}

\item{diagonal}{for \code{diagonal=FALSE} (the default) CAR scores are computed;
                otherwise with \code{diagonal=TRUE} marginal correlations.
	     }

\item{verbose}{if \code{verbose=TRUE} then the shrinkage intensity 
                used in estimating the shrinkage correlation matrix is reported.
	     }


}
\description{
  \code{carscore} estimates the CAR score vector, either using the standard
  empirical estimator of the correlation matrix, or a shrinkage estimator.
}

\value{
  \code{carscore} returns
  a vector containing the CAR scores (or marginal correlations).
}

\details{
   The CAR scores are the correlations between the response and the 
   Mahalanobis-decorrelated predictors.  In Zuber and Strimmer (2010) it
   is argued that squared CAR scores are a natural measure for variable 
   importance --- see \url{http://arxiv.org/abs/1007.5516} for details.
}

\author{
  Verena Zuber and Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link{care}}.}


\references{

Zuber, V., and K. Strimmer. 2010.  Variable importance and model selection by decorrelation.
\url{http://arxiv.org/abs/1007.5516}.
}


\examples{
# load care library
library("care")

######

# empirical CAR scores for diabetes data
data(efron2004)
car = carscore(efron2004$x, efron2004$y)
car

# R^2
sum(car^2)

# variables in order of squared CAR scores
xnames = colnames(efron2004$x)
xnames[order(car^2, decreasing=TRUE)]

# compare with standard ordering given by t-scores / partial correlations 
fit = fitlm(efron2004$x, efron2004$y)
xnames[order(fit$pcor^2, decreasing=TRUE)]

######

# shrinkage CAR scores for Lu et al. (2004) data
data(lu2004)
dim(lu2004$x)    # 30 403

# compute shrinkage car score
car = carscore(lu2004$x, lu2004$y, estimator="shrinkage")

# R^2
sum(car^2) 

# most important genes
order(car^2, decreasing=TRUE)[1:10]

# compare with empirical marginal correlations
marg = carscore(lu2004$x, lu2004$y, diagonal=TRUE)
order(abs(marg), decreasing=TRUE)[1:10]

}
\keyword{multivariate}
