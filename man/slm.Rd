\name{slm}
\alias{slm}
\alias{predict.slm}


\title{Shrinkage Estimation of Regression Coefficients}
\description{
  \code{slm} fits a linear model and computes
 (standardized) regression coefficients by plugin of shrinkage estimates of correlations and variances.
 
}

\usage{
slm(Xtrain, Ytrain, lambda, lambda.var, diagonal=FALSE, verbose=TRUE)
\method{predict}{slm}(object, Xtest, verbose=TRUE, ...)
}

\arguments{
   \item{Xtrain}{Matrix of predictors (columns correspond to variables).}
  
  \item{Ytrain}{Univariate continous response variable.}
  
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

  \code{diagonal}: Whether or not the correlations among predictors were set to zero.

  \code{std.coefficients}: The standardized regression coefficients, i.e. the regression coefficients
           computed from centered and standardized input data.  Thus, by construction the intercept is zero.
           Furthermore, for \code{diagonal=TRUE} the standardized regression coefficient for each predictor is
           identical to the respective marginal correlation.

  \code{coefficients}: Regression coefficients.

  \code{R2}: For \code{diagonal=TRUE} this is the multiple correlation coefficient 
             between the response and the predictor, or the proportion of explained variance, with range
             from 0 to 1.
             For \code{diagonal=TRUE} this equals the sum of squared marginal
              correlations.  Note that this sum may be larger than 1!
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


Sch\"afer, J., and K. Strimmer. 2005.  A shrinkage approach to
 large-scale covariance estimation and implications for functional
 genomics.  Statist. Appl. Genet. Mol. Biol. 4: 32.  
 (\url{http://www.bepress.com/sagmb/vol4/iss1/art32/})
}
\seealso{\code{\link{carscore}}, \code{\link{slm.models}}}

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
b = fit$coefficients[-1]
cbind(b, pcor, t, t.pval)

# compare results with those from lm function
lm.out = lm(y ~ x)
summary(lm.out)

# prediction of fitted values at the position of the training data
predict(fit, x)
lm.out$fitted.values


## example with small number of samples and large dimension
## (using shrinkage estimates of regression coefficients)

data(lu2004)
dim(lu2004$x)    # 30 403

fit = slm(lu2004$x, lu2004$y)
fit
fit$regularization

predict(fit, lu2004$x)-lu2004$y

}
\keyword{multivariate}
