\name{ProfileScore}
\alias{ProfileScore}
\title{Profile Parsimony score}
\description{\dots}
\usage{
ProfileScore (tree, data)
}
\arguments{
  \item{tree}{An object of class \code{phyDat} denoting the topology to begin the search from;}
  \item{data}{the output of \code{\link{PrepareDataFitch}}; alternatively the function will 
              convert (time-consumingly) data of class \code{phyDat}.}
}
\details{
    \dots
}
\value{
  The Profile Parsimony score, expressed as a negative value so that smaller values are better.
}
\references{
  \itemize{
    \item Faith, D. P. & Trueman, J. W. H. (2001). \cite{Towards an inclusive philosophy for phylogenetic
    inference.} Systematic Biology 50:3, 331-350, doi:
    \href{http://dx.doi.org/10.1080/10635150118627}{10.1080/10635150118627}
  }
}

\author{
Martin R. Smith

Adapted from \code{\link{pratchet}} in the \pkg{phangorn} package, which does not preserve the position of the root.
}
\seealso{
\itemize{
\item \code{\link{Pratchet}}
\item \code{\link{TreeSearch}}
}
}
\examples{
data(referenceTree)
data(congreveLamsdellMatrices)
dataset <- PrepareDataFitch (congreveLamsdellMatrices[[42]])
ProfileScore(referenceTree, dataset)
}
\keyword{ tree }