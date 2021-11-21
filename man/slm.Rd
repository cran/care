\name{slm}
\alias{slm}
\alias{predict.slm}
\alias{make.predlist}


\title{Shrinkage Estimation of Regression Coefficients}
\description{
  \code{slm} fits a linear model and computes
 (standardized) regression coefficients by plugin of shrinkage estimates of correlations and variances.
 Using the argument \code{predlist} several models can be fitted on the same data set.
 
 \code{make.predlist} constructs a \code{predlist} argument for use with \code{slm}.

}

\usage{
slm(Xtrain, Ytrain, predlist, lambda, lambda.var, diagonal=FALSE, verbose=TRUE)
\method{predict}{slm}(object, Xtest, verbose=TRUE, ...)
make.predlist(ordering, numpred, name="SIZE")
}

\arguments{
   \item{Xtrain}{Matrix of predictors (columns correspond to variables).}
  
  \item{Ytrain}{Univariate continous response variable.}

 \item{predlist}{A list specifying the predictors to be included when fitting the linear regression.  Each entry in the list is a vector containing the indices of variables used per model. If left unspecified single full-sized model using all variables in Xtrain is assumed. For a given ordering of covariables a suitable \code{predlist} can be generated  using the helper function \code{\link{make.predlist}} - see examples below.
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

  \item{object}{An \code{slm} fit object obtained from the function \code{slm}.}
  
  \item{Xtest}{A matrix containing the test data set. Note that 
                the rows correspond to observations and the columns
                to variables.}
    \item{...}{Additional arguments for generic predict.}

 \item{ordering}{The ordering of the predictors (most important predictors are first).
	     }


 \item{numpred}{The number of included predictors (may be a scalar or a vector). The predictors
               are included in the order specified by \code{ordering}.
	     }

 \item{name}{The name assigned to each model is \code{name} plus "." and the number of included predictors.
	     }
  
}
\details{
 
 The regression coefficients are obtained by estimating the joint
 joint covariance matrix of the response and the predictors,
 and subsequently computing the the regression coefficients by
 inversion of this matrix - see  Opgen-Rhein and Strimmer (2007).
 As estimators for the covariance matrix either the standard empirical estimator
 or a Stein-type shrinkage estimator is employed.  The use of the empirical covariance
 leads to the OLS estimates of the regression coefficients, whereas otherwise
 shrinkage estimates are obtained. 
}
\value{
  \code{slm} returns a list with the following components:

  \code{regularization}: The shrinkage intensities used for estimating correlations and variances.

  \code{std.coefficients}: The standardized regression coefficients, i.e. the regression coefficients
           computed from centered and standardized input data.  Thus, by construction the intercept is zero.
           Furthermore, for \code{diagonal=TRUE} the standardized regression coefficient for each predictor is
           identical to the respective marginal correlation.

  \code{coefficients}: Regression coefficients.

  \code{numpred}: The number of predictors used in each investigated model.


  \code{R2}: For \code{diagonal=TRUE} this is the multiple correlation coefficient 
             between the response and the predictor, or the proportion of explained variance, with range
             from 0 to 1.
             For \code{diagonal=TRUE} this equals the sum of squared marginal
              correlations.  Note that this sum may be larger than 1!

  \code{sd.resid}: The residual unexplained error.

  \code{predict.slm} returns the means predicted for each sample and model as well as the corresponding
         predictive standard deviations (attached as attribute "sd").


}
\author{
  Korbinian Strimmer (\url{https://strimmerlab.github.io}).
  
}
\references{

 Opgen-Rhein, R., and K. Strimmer. 2007.
 From correlation to causation networks: a simple approximate learning 
 algorithm and its   application to high-dimensional plant gene expression data.
 BMC Syst. Biol. 1: 37.  <DOI:10.1186/1752-0509-1-37>

Sch\"afer, J., and K. Strimmer. 2005.  A shrinkage approach to
 large-scale covariance estimation and implications for functional
 genomics.  Statist. Appl. Genet. Mol. Biol. 4: 32.  <DOI:10.2202/1544-6115.1175>
}
\seealso{\code{\link{carscore}} }

\examples{
# load care library
library("care")

## example with large number of samples and small dimension
## (using empirical estimates of regression coefficients)

# diabetes data
data(efron2004)
x = efron2004$x
y = efron2004$y
n = dim(x)[1]
d = dim(x)[2]
xnames = colnames(x)

# empirical regression coefficients
fit = slm(x, y, lambda=0, lambda.var=0)
fit
# note that in this example the regression coefficients
# and the standardized regression coefficients are identical
# as the input data have been standardized to mean zero and variance one

# compute corresponding t scores / partial correlations
df = n-d-1
pcor = pcor.shrink(cbind(y,x), lambda=0)[-1,1] 
t = pcor * sqrt(df/(1-pcor^2))
t.pval = 2 - 2 * pt(abs(t), df)
b = fit$coefficients[1,-1]
cbind(b, pcor, t, t.pval)

# compare results with those from lm function
lm.out = lm(y ~ x)
summary(lm.out)

# prediction of fitted values at the position of the training data
lm.out$fitted.values
mu.hat = predict(fit, x) # precticted means
mu.hat
attr(mu.hat, "sd") # predictive error
sd(y-mu.hat)


# ordering of the variables using squared empirical CAR score
car = carscore(x, y, lambda=0)
ocar = order(car^2, decreasing=TRUE)
xnames[ocar]

# CAR regression models with 5, 7, 9 included predictors
car.predlist = make.predlist(ocar, numpred = c(5,7,9), name="CAR")
car.predlist
slm(x, y, car.predlist, lambda=0, lambda.var=0)


# plot regression coefficients for all possible CAR models

p=ncol(x)
car.predlist = make.predlist(ocar, numpred = 1:p, name="CAR")
cm = slm(x, y, car.predlist, lambda=0, lambda.var=0)
bmat = cm$coefficients[,-1]
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


## example with small number of samples and large dimension
## (using shrinkage estimates of regression coefficients)

data(lu2004)
dim(lu2004$x)    # 30 403

fit = slm(lu2004$x, lu2004$y)
fit

}
\keyword{multivariate}
