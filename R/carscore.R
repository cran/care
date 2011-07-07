### carscore.R  (2011-07-03)
###
###    Estimate CAR scores and marginal correlations
###
### Copyright 2010-2011 Verena Zuber and Korbinian Strimmer
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




# estimate car scores 
# (and marginal correlations if diagonal=TRUE)

carscore = function(Xtrain, Ytrain, diagonal=FALSE, shrink=TRUE, verbose=TRUE)
{
  n = dim(Xtrain)[1]
  p = dim(Xtrain)[2]

  #####################################

  omega = cor(Xtrain, Ytrain)  # marginal correlations
  if (shrink) # shrinkage estimator
  {
    # regularize the joint correlation matrix  Ytrain and Xtrain combined
    if(verbose) cat("Estimating optimal shrinkage intensity lambda (correlation matrix): ")
    lambda = pvt.corlambda(scale(cbind(Ytrain,Xtrain)), rep(1/n, n), 0)
    if(verbose) cat(round(lambda, 4), "\n")

    omega = (1-lambda)*omega # shrink marginal correlations
  }
  else # empirical estimator
  {
    lambda = 0
  }

  if (diagonal==FALSE)
  {
      # car score
      omega = crossprod.powcor.shrink(Xtrain, omega, alpha=-1/2, lambda=lambda, verbose=FALSE)
  }

  omega = as.vector(omega)
  names(omega) = colnames(Xtrain)
  if(shrink)
  {
    class(omega) = "shrinkage"
    attr(omega, "lambda") = lambda
    attr(omega, "lambda.estimated") = TRUE
  }

  return( omega )
}

