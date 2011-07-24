\name{slm.models}
\alias{slm.models}
\alias{make.predlist}

\title{Compare Different Regression Models}
\usage{
slm.models(Xtrain, Ytrain, predlist, lambda, lambda.var, diagonal=FALSE, verbose=TRUE)
make.predlist(ordering, numpred, name="SIZE")
}
\arguments{
  \item{Xtrain}{Matrix of predictors (columns correspond to variables).}
  
  \item{Ytrain}{Univariate response variable.}

 \item{predlist}{A list specifying the predictors included in each
                 model. If the predictors are ordered such a list can be generated
                 using the function \code{\link{make.predlist}}.
	     }

  \item{lambda}{The correlation shrinkage intensity (range 0-1).  
          If not specified (the default) it is estimated using an
          analytic formula from Sch\"afer and Strimmer (2005).  For \code{lambda=0}
      the empirical correlations are used.}

  \item{lambda.var}{The variance shrinkage intensity (range 0-1).  If
           not specified (the default) it is estimated
          using an analytic formula from Opgen-Rhein and Strimmer
          (2007).  For \code{lambda.var=0} the empirical
          variances are used.}

  \item{diagonal}{If \code{diagonal=FALSE} (the default) then the correlation
                  among predictor veriables assumed to be non-zero and is estimated
                  from data. If \code{diagonal=TRUE} 
                  then it is assumed that the correlation among predictors vanishes
                  and is set to zero.}

  \item{verbose}{If \code{verbose=TRUE} then the estimated shrinkage intensities are reported.
	     }


  \item{ordering}{The ordering of the predictors (most important predictors are first).
	     }


 \item{numpred}{The number of included predictors (may be a scalar or a vector). The predictors
               are included in the order specified by \code{ordering}.
	     }

 \item{name}{The name assigned to each model is \code{name} plus "." and the number of included predictors.
	     }
}
\description{
  \code{slm.models} estimates regression coefficients and R2 for all the specified models.
  This is useful for comparing different regression models (see examples below).

  \code{make.predlist} constructs a \code{predlist} argument for use with \code{slm.models}.
}

\value{
  \code{slm.models} returns a list with the following components:

  \code{coefficients}: Regression coefficients for all investigated models.

  \code{numpred}: The number of predictors used in each investigated model.

  \code{R2}: The value of R2 for all investigated models.

}

\author{
  Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link{slm}}, \code{\link{carscore}}.}


\references{

Zuber, V., and K. Strimmer. 2011.  High-dimensional regression and 
variable selection using CAR scores.  Statist. Appl. Genet. Mol. Biol. 10: 34.
\url{http://www.bepress.com/sagmb/vol10/iss1/art34/}

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
car = carscore(x, y, lambda=0)
ocar = order(car^2, decreasing=TRUE)
xnames[ocar]

# CAR regression models with 5, 7, 9 included predictors
car.predlist = make.predlist(ocar, numpred = c(5,7,9), name="CAR")
car.predlist
slm.models(x, y, car.predlist, lambda=0, lambda.var=0)


#####

# plot regression coefficients for all possible CAR models

p=ncol(x)
car.predlist = make.predlist(ocar, numpred = 1:p, name="CAR")
cm = slm.models(x, y, car.predlist, lambda=0, lambda.var=0)
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
