### predict.slm.R  (2011-07-03)
###
###    Prediction from linear model
###
### Copyright 2011 Korbinian Strimmer
###
###
### This file is part of the `sda' library for R and related languages.
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




predict.slm = function(object, Xtest, verbose=TRUE, ...)
{
  if ( missing(object) ) {
    stop("An slm object must be supplied.")
  }

  if ( missing(Xtest) ) {
    stop("A test data set must be supplied.")
  }
  
  if (!is.matrix(Xtest)) stop("Test data must be given as matrix!")
  ntest = nrow(Xtest) # number of test samples
  npredTest = ncol(Xtest) # number of predictor variables in test data set
  npred =  length(object$coefficients)-1 # number of predictor variables 

  if (npred != npredTest)
    stop("Different number of predictors in slm object (", 
         npred, ") and in test data (", 
         npredTest, ")", sep="")

  if (verbose) cat("Prediction uses", npred, "variables.\n")

  b = matrix(object$coefficients[-1]) 
  b0 = object$coefficients[1]

  yhat =  b0 + Xtest %*% b 
 
  return( yhat )
}

