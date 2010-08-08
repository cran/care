### care.R  (2010-08-06)
###
###    CAR regression models
###
### Copyright 2010 Korbinian Strimmer
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


# CAR regression models

care = function(x, y, numpred, estimator=c("empirical", "shrinkage"), verbose=TRUE)
{
  p = dim(x)[2]
  m = length(numpred)

  coeff = matrix(0, nrow=m, ncol=p+1) 
  rownames(coeff) = paste("(numpred=", numpred, ")", sep="")
  xnames = colnames(x)
  if ( is.null(xnames) ) 
    xnames = paste("X", 1:p, sep="")
  colnames(coeff) = c("(Intercept)", xnames)
  
  
  # estimate CAR score
  if(verbose) cat("Determine CAR ordering of the variables\n")
  car = carscore(x,y, estimator=estimator, verbose=verbose)
  # order variables by squared CAR scores
  ocar = order(car^2, decreasing=TRUE)
  if(verbose) cat("\nDetermine regression coefficients\n")
  

  for (i in 1:m)
  {
    # keep the best variables
    idx = ocar[1:(numpred[i])]

    fit = fitlm(x[,idx, drop=FALSE], y, estimator=estimator, verbose=verbose)
    coeff[i, 1+idx] = fit$b
    coeff[i, 1] = fit$a
  }

  return(coeff)
}



