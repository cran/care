\name{lu2004}
\alias{lu2004}

\title{Brain Aging Study of Lu et al. (2004)}


\usage{
data(lu2004)
}

\description{
Gene expression data (403 genes for 30 samples) from the
microarray study of Lu et al. (2004).  
}

\details{ This data set contains measurements of the gene
expression of 403 genes from 30 human brain samples.  In addition,
the age of each patient is provided.  
}

\format{
  \code{lu2004$x} is a 30 x 403 matrix containing the log expression levels.
   The rows contain the samples and the columns the genes.

   \code{lu2004$y} is the age of for each sample.
}

\source{
The orginal data are available from \url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE1572} and are described in Lu et al. (2004).
The selected 403 genes result from prescreening and preprocessing as
described in Zuber and Strimmer (2010). 
}

\references{
Lu, T., et al. 2004. Gene regulation and DNA damage in the ageing human brain.
Nature 429:883--891.

Zuber, V., and K. Strimmer. 2010.  Variable importance and model selection by decorrelation.
\url{http://arxiv.org/abs/1007.5516}.
}


\examples{
# load care library
library("care")

# load Lu et al. (2004) data set
data(lu2004)
dim(lu2004$x) # 30 403
hist(lu2004$x)
length(lu2004$y) # 30
lu2004$y # age
}


\keyword{datasets}
