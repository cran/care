\name{efron2004}
\alias{efron2004}

\title{Diabetes Data from Efron et al. (2004)}


\usage{
data(efron2004)
}

\description{
Diabetes data (10 variables, 442 measurements) as used in 
the study of Efron et al. (2004).  The data is standardized such that
the means of all variables are zero, and all variances are equal to one.
}

\format{
  \code{efron2004$x} is a 422 x 10 matrix containing the measurements
    of the explanatory variables (age, sex, body mass, etc.).
   The rows contain the samples and the columns the variables.

   \code{efron2004$y} contains the response.
}

\source{
The orginal data are available in the lars R package, 
see \url{https://cran.r-project.org/package=lars}. Note that this uses a slightly different standardization.
}

\references{
Efron, B., et al. 2004. Least angle regression (with discussion). Ann. Statist. 32:407--499. <DOI:10.1214/009053604000000067>
}


\examples{
# load care library
library("care")

# load Efron et al. (2004) diabetes data set
data(efron2004)
dim(efron2004$x) # 442 10
colnames(efron2004$x)
length(efron2004$y) # 442
}


\keyword{datasets}
