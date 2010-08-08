\name{fitlm}
\alias{fitlm}


\title{Estimation of Regression Coefficients}
\description{
  \code{fitlm} fits regression coefficients by plugin of shrinkage 
    or empirical estimates of correlations and variances.
}

\usage{
fitlm(x, y, estimator=c("empirical", "shrinkage"), verbose=TRUE)
}

\arguments{
   \item{x}{matrix of predictors (columns correspond to variables).}
  
  \item{y}{univariate response variable.}

  \item{estimator}{either \code{"empirical"} (the default) or \code{"shrinkage"}.}

  \item{verbose}{if \code{verbose=TRUE} then the shrinkage intensity 
                used in estimating the shrinkage correlation matrix is reported.
	     }

  
}
\details{
 
 The regression coefficients are fitted by estimating the joint covariance
 matrix of the response and the predictors, which is then inverted and
 used for computing regression coefficients - see  Opgen-Rhein and Strimmer (2007).
 As estimators for the covariance matrix either the standard empirical estimator
 or a Stein-type shrinkage estimator is employed.  The use of the empirical covariance
 leads to the OLS estimates of the regression coefficients, whereas otherwise
 shrinkage estimates are obtained. 
 
 The function also returns the partial 
 partial correlations between each predictor and the response.
 These provide the standard ordering
 of variables (identical to the ordering of p-values or t-scores computed by    
 \code{lm}, see example below) but different from the ranking provided
 by CAR scores (cf.  \code{\link{carscore}}).

}
\value{
  \code{fitlm} returns a list with the estimated partial correlations, the
   regression coefficients, and the intercept. 
}
\author{
  Korbinian Strimmer (\url{http://strimmerlab.org}).
  
}
\references{

 Opgen-Rhein, R., and K. Strimmer. 2007.
 From correlation to causation networks: a simple approximate learning 
 algorithm and its   application to high-dimensional plant gene expression data.
 BMC Syst. Biol. 1: 37. 
 (\url{http://www.biomedcentral.com/1752-0509/1/37/abstract})

}
\seealso{\code{\link{carscore}}}

\examples{
# load care library
library("care")

# diabetes data
data(efron2004)
x = efron2004$x
y = efron2004$y
n = dim(x)[1]
p = dim(x)[2]

# empirical regression coefficients
fit = fitlm(x, y, estimator="empirical")
fit

# t scores corresponding to partial correlations
df = n-p-1
r = fit$pcor
t = r * sqrt(df/(1-r^2))
t.pval = 2 - 2 * pt(abs(t), df)

cbind(fit$b, t, t.pval)

# compare results with lm
summary(lm(y ~ x))

}
\keyword{multivariate}
