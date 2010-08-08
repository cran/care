\name{care}
\alias{care}

\title{Compute CAR Regression Models}
\usage{
care(x, y, numpred, estimator=c("empirical", "shrinkage"), verbose=TRUE)
}
\arguments{
  \item{x}{matrix of predictors (columns correspond to variables).}
  
  \item{y}{univariate response variable.}

\item{numpred}{number of included predictors (may be a scalar or a vector).
	     }


  \item{estimator}{either \code{"empirical"} (the default) or \code{"shrinkage"}.}


\item{verbose}{if \code{verbose=TRUE} then the shrinkage intensity 
                used in estimating the shrinkage correlation matrix is reported.
	     }


}
\description{
  \code{care} estimates CAR regression models, i.e. it first computes the 
  CAR scores. Subsequently, the variables are included in the model in the order given
  by the squared CAR scores.  See Zuber and Strimmer (2010) for details.
}

\value{
  \code{care} returns a matrix of regression coefficients.
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
\seealso{\code{\link{carscore}}, \code{\link{fitlm}}.}


\references{

Zuber, V., and K. Strimmer. 2010.  Variable importance and model selection by decorrelation.
\url{http://arxiv.org/abs/1007.5516}.
}


\examples{
library("care")

# diabetes data
data(efron2004)
x = efron2004$x
y = efron2004$y
xnames = colnames(x)

#####

# ordering of the variables using squared CAR score
car = carscore(x, y)
xnames[order(car^2, decreasing=TRUE)]

# CAR regression models with 5, 7, 9 included predictors
care(x, y, numpred=c(5,7,9))

#####

# plot regression coefficients for all possible CAR models
p = 10
bmat= care(x, y, 1:p)[,-1]
bmat

plot(1:p, bmat[,1], type="l", 
  ylab="estimated regression coefficients", 
  xlab="number of included predictors", 
  main="CAR Regression Models for Diabetes Data", 
  xlim=c(1,p+1), ylim=c(min(bmat), max(bmat)))

for (i in 2:p) lines(1:p, bmat[,i], col=i, lty=i)
for (i in 1:p) points(1:p, bmat[,i], col=i)
for (i in 1:p) text(p+0.5, bmat[p,i], xnames[i])


}
\keyword{multivariate}
