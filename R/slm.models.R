### slm.models.R  (2011-07-24)
###
###    Compare different regression models
###
### Copyright 2010-11 Korbinian Strimmer
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


# generic function

slm.models = function(Xtrain, Ytrain, predlist, lambda, lambda.var, diagonal=FALSE, verbose=TRUE)
{
  p = dim(Xtrain)[2]
  m = length(predlist)
  numpred = sapply(predlist, length)
  R2 = numeric(m)
  coeff = matrix(0, nrow=m, ncol=p+1)
  
  modelnames = names(predlist)
  if ( is.null(modelnames) ) 
    modelnames = paste("SIZE.", numpred, sep="")
  xnames = colnames(Xtrain)
  if ( is.null(xnames) ) 
    xnames = paste("X", 1:p, sep="")

  rownames(coeff) = modelnames  
  colnames(coeff) = c("(Intercept)", xnames)
  names(R2) = modelnames
  names(numpred) = modelnames
    
  for (i in 1:m)
  {
   if(verbose) cat("Determine regression coefficients for", modelnames[i], "model\n")

    # keep the best variables
    idx = predlist[[i]]

    fit = slm(Xtrain[,idx, drop=FALSE], Ytrain, diagonal=diagonal, lambda=lambda, lambda.var=lambda.var, verbose=verbose)

    coeff[i, 1+idx] = fit$coefficients[-1] # predictors
    coeff[i, 1] = fit$coefficients[1] # intercept
    R2[i] = fit$R2
  }

  return(list(coefficients=coeff, numpred=numpred, R2=R2))
}


make.predlist = function(ordering, numpred, name="SIZE")
{
  predlist = vector("list", length(numpred))
  names(predlist) =  paste(name, ".", numpred, sep="")
  for (i in 1:length(numpred))
    predlist[[i]] = ordering[1:numpred[i]]

  return(predlist)
}

