
### fitlm.R  (2010-08-06)
###
###    Fit regression coefficients by plugin of (shrinkage or empirical) 
###    estimates of correlations and variances
###
### Copyright 2006-2010 Korbinian Strimmer
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



fitlm = function(x, y, estimator=c("empirical", "shrinkage"), verbose=TRUE)
{
  estimator = match.arg(estimator)
  
  n = dim(x)[1]
  p = dim(x)[2]
  yx = cbind(y,x)
  mu = apply(yx, 2, mean)
  
  if ( n <= p+1  & estimator=="empirical")
  {
    warning("Singular empirical estimator, using shrinkage estimator instead.\n") 
    estimator="shrinkage"
  }

  if (estimator=="empirical")
  {
    prec = solve( cov(yx) )
  }
 
  if (estimator=="shrinkage")
  {
    prec = invcov.shrink(yx, verbose=verbose)
  }

  pcor =-prec[1,-1]/sqrt(prec[1,1]*diag(prec)[-1])  # partial correlations
  b =-prec[1,-1]/prec[1,1]                          # regression coefficients
  a = mu[1] - sum(b * mu[-1])                       # intercept

  names(a) = c("(Intercept)")
  return( list(pcor=pcor, b=b, a=a) )
}


