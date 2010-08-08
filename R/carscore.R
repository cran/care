### carescore.R  (2010-08-06)
###
###    Estimate CAR scores and marginal correlations
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




# estimate car scores 
# (and marginal correlations if diagonal=TRUE)

carscore = function(x, y, estimator=c("empirical", "shrinkage"), 
               diagonal=FALSE, verbose=TRUE)
{
  estimator = match.arg(estimator)
  
  n = dim(x)[1]
  p = dim(x)[2]
  
  if ( n <= p+1 & estimator=="empirical" & diagonal == FALSE)
  {
    warning("Singular empirical estimator, using shrinkage estimator instead.\n") 
    estimator="shrinkage"
  }

  #####################################

  if (estimator=="empirical")
  {
    Rxy = cor(x, y)
    if (diagonal==FALSE)
    {
      Rxx.invroot = mpower(cor(x), -1/2)
    }
  }

  if (estimator=="shrinkage")
  {
    # regularize the joint correlation matrix  y and x combined
    w = rep(1/n, n)
    if(verbose) cat("Estimating optimal shrinkage intensity lambda (correlation matrix): ")
    lambda = pvt.corlambda(scale(cbind(x,y)), w, 0)
    if(verbose) cat(round(lambda, 4), "\n")

    Rxy = (1-lambda)* cor(x, y)
    if (diagonal==FALSE)
    {
      Rxx.invroot = powcor.shrink(x, alpha=-1/2, 
                        lambda=lambda, w=w, verbose=FALSE)
    }
  }

  if (diagonal==FALSE)
    omega = as.vector( Rxx.invroot %*% Rxy )   # car score
  else
    omega = as.vector( Rxy )                   # marginal correlations

  names(omega) = colnames(x)

  return( omega )
}

