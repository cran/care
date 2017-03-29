\name{carscore}
\alias{carscore}

\title{Estimate CAR Scores and Marginal Correlations}
\usage{
carscore(Xtrain, Ytrain, lambda, diagonal=FALSE, verbose=TRUE)
}
\arguments{
  \item{Xtrain}{Matrix of predictors (columns correspond to variables).}
  
  \item{Ytrain}{Univariate response variable.}

  \item{lambda}{The correlation shrinkage intensity (range 0-1).  
          If not specified (the default) it is estimated using an
          analytic formula from Sch\"afer and Strimmer (2005).  For \code{lambda=0}
      the empirical correlations are used.}

\item{diagonal}{For \code{diagonal=FALSE} (the default) CAR scores are computed;
                otherwise with \code{diagonal=TRUE} marginal correlations.   
	     }

\item{verbose}{If \code{verbose=TRUE} then the shrinkage intensity 
                used in estimating the shrinkage correlation matrix is reported.
	     }


}
\description{
  \code{carscore} estimates the vector of CAR scores, either using the standard
  empirical estimator of the correlation matrix, or a shrinkage estimator.
}

\value{
  \code{carscore} returns
  a vector containing the CAR scores (or marginal correlations for \code{diagonal=TRUE}).
}

\details{
   The CAR scores are the correlations between the response and the 
   Mahalanobis-decorrelated predictors.  CAR score is an abbreviation
   for Correlation-Adjusted (marginal) coRelation, where the first 
   correlation matrix refers dependencies among predictors.

    In Zuber and Strimmer (2011) it
   is argued that squared CAR scores are a natural measure for variable 
   importance and it is shown that variable selection based on CAR scores
   is highly efficient compared to competing approaches such as elastic net
   lasso, or boosting. 

   If the response is binary (or descrete) the corresponding quantity
   are CAT scores (see \code{\link[sda]{catscore}}).
}

\author{
  Verena Zuber and Korbinian Strimmer (\url{http://strimmerlab.org}).
}
\seealso{\code{\link[sda]{catscore}}.}


\references{

Zuber, V., and K. Strimmer. 2011.  High-dimensional regression and 
variable selection using CAR scores.  Statist. Appl. Genet. Mol. Biol. 10: 34.
<DOI:10.2202/1544-6115.1730>
}


\examples{
# load care library
library("care")

######

# empirical CAR scores for diabetes data
data(efron2004)
xnames = colnames(efron2004$x)
n = dim(efron2004$x)[1]

car = carscore(efron2004$x, efron2004$y, lambda=0)
car

# compare orderings

# variables ordered by squared CAR scores
xnames[order(car^2, decreasing=TRUE)]
# "bmi" "s5"  "bp"  "s3"  "s4"  "s6"  "sex" "age" "s2"  "s1" 

# compare with ordering by t-scores / partial correlations 
pcor = pcor.shrink(cbind(efron2004$y,efron2004$x), lambda=0, verbose=FALSE)[-1,1]
xnames[order(pcor^2, decreasing=TRUE)]
# "bmi" "bp"  "s5"  "sex" "s1"  "s2"  "s4"  "s6"  "s3"  "age"

# compare with ordering by marginal correlations 
mcor = cor(efron2004$y,efron2004$x)
#mcor = carscore(efron2004$x, efron2004$y, diagonal=TRUE, lambda=0)
xnames[order(mcor^2, decreasing=TRUE)]
# "bmi" "s5"  "bp"  "s4"  "s3"  "s6"  "s1"  "age" "s2"  "sex"

# decomposition of R^2
sum(car^2)
slm(efron2004$x, efron2004$y, lambda=0, lambda.var=0)$R2

# pvalues for empirical CAR scores
pval = 1-pbeta(car^2, shape1=1/2, shape2=(n-2)/2)
pval <= 0.05

######

# shrinkage CAR scores for Lu et al. (2004) data
data(lu2004)
dim(lu2004$x)    # 30 403

# compute shrinkage car scores
car = carscore(lu2004$x, lu2004$y)

# most important genes
order(car^2, decreasing=TRUE)[1:10]

# compare with empirical marginal correlations
mcor = cor(lu2004$y, lu2004$x)
order(mcor^2, decreasing=TRUE)[1:10]

# decomposition of R^2
sum(car^2) 
slm(lu2004$x, lu2004$y)$R2

}
\keyword{multivariate}
