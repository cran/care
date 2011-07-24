
### slm.R  (2011-07-23)
###
###    Fit regression coefficients by plugin of (shrinkage or empirical) 
###    estimates of correlations and variances
###
### Copyright 2006-2011 Korbinian Strimmer
###
###
### This file is part of the `care' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA



slm = function(Xtrain, Ytrain, lambda, lambda.var, diagonal=FALSE, verbose=TRUE)
{
  n = dim(Xtrain)[1]
  p = dim(Xtrain)[2]
  yx = cbind(Ytrain,Xtrain)
  mu = apply(yx, 2, mean)

  # regularize the joint correlation matrix  y and x combined
  if(missing(lambda))
  {
    if(verbose) cat("Estimating optimal shrinkage intensity lambda (correlation matrix): ")
    lambda = pvt.corlambda(scale(yx), rep(1/n, n), 0)
    if(verbose) cat(round(lambda, 4), "\n")
  }
  else
  {
     if(verbose) cat("Specified shrinkage intensity lambda (correlation matrix): ", round(lambda, 4), "\n")
  }

  mcor = (1-lambda)*cor(Xtrain, Ytrain) # marginal correlations

  v = var.shrink(yx, lambda.var=lambda.var, verbose=verbose)
  lambda.var = attr(v, "lambda.var")
  
  regularization = c(lambda, lambda.var)
  names(regularization) = c("lambda", "lambda.var")
  sdy = sqrt(v[1])
  sdx = sqrt(v[-1])

  if (diagonal)
  {
    bstd = mcor
  }
  else
  {
    bstd = crossprod.powcor.shrink(Xtrain, mcor, alpha=-1, lambda=lambda, verbose=FALSE)
  }

  R2 = as.vector(crossprod(mcor, bstd))   # proportion of explained variance  (diagonal=FALSE)
                                          # or sum of squared marginal correlations (diagonal=TRUE)

  b = bstd*sdy/sdx                   # regression coefficients
  a = mu[1] - sum(b * mu[-1])        # intercept

  coefficients = c(a, b)
  xnames = colnames(Xtrain)
  if ( is.null(xnames) ) 
    xnames = paste("X", 1:p, sep="")
  names(coefficients) = c("(Intercept)", xnames)

  std.coefficients = c(0, bstd)
  names(std.coefficients) = c("(Intercept)", xnames)

  res = list(regularization=regularization, diagonal=diagonal, 
          std.coefficients=std.coefficients,
          coefficients=coefficients, R2=R2) 
  class(res) = "slm"

  return( res )
}

