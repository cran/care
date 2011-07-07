\name{car.models}
\alias{car.models}

\title{Compute CAR Regression Models}
\usage{
car.models(x, y, numpred, shrink=TRUE, verbose=TRUE)
}
\arguments{
  \item{x}{matrix of predictors (columns correspond to variables).}
  
  \item{y}{univariate response variable.}

\item{numpred}{number of included predictors (may be a scalar or a vector). The predictors
               are included in the order of squared CAR scores.
	     }


  \item{shrink}{Use a shrinkage estimator for estimation of correlations and variances
                (default). Alternatively, employ the standard empirical estimates.}


\item{verbose}{if \code{verbose=TRUE} then the shrinkage intensity 
                used in estimating the shrinkage correlation matrix is reported.
	     }


}
\description{
  \code{car.models} estimates CAR regression models. First, the CAR scores
  using the whole data set are computed.
   Subsequently, the variables are included in the model in the order given
  by the squared CAR scores, and regression coefficients are computed. 
}

\value{
  \code{car.models} returns a matrix of regression coefficients and the corresponding estimated R2s.
}

\author{
  Verena Zuber and Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link{carscore}}, \code{\link{slm}}.}


\references{

Zuber, V., and K. Strimmer. 2011.  High-dimensional regression and 
variable selection using CAR scores.  Statist. Appl. Genet. Mol. Biol. 10: in press. 
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

# ordering of the variables using squared empirical CAR score
car = carscore(x, y, shrink=FALSE)
xnames[order(car^2, decreasing=TRUE)]

# CAR regression models with 5, 7, 9 included predictors
car.models(x, y, numpred=c(5,7,9))

#####

# plot regression coefficients for all possible CAR models
p = 10
cm = car.models(x, y, numpred=1:p, shrink=FALSE)
bmat= cm$coefficients[,-1]
bmat

par(mfrow=c(2,1))

plot(1:p, bmat[,1], type="l", 
  ylab="estimated regression coefficients", 
  xlab="number of included predictors", 
  main="CAR Regression Models for Diabetes Data", 
  xlim=c(1,p+1), ylim=c(min(bmat), max(bmat)))

for (i in 2:p) lines(1:p, bmat[,i], col=i, lty=i)
for (i in 1:p) points(1:p, bmat[,i], col=i)
for (i in 1:p) text(p+0.5, bmat[p,i], xnames[i])

plot(1:p, cm$R2, type="l", 
  ylab="estimated R2",
  xlab="number of included predictors",
  main="Proportion of Explained Variance",
  ylim=c(0,0.6))
R2max = max(cm$R2)
lines(c(1,p), c(R2max, R2max), col=2)

par(mfrow=c(1,1))

}
\keyword{multivariate}
