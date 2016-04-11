\name{referenceTree}
\alias{referenceTree}
\docType{data}
\title{Reference tree from Congreve and Lamsdell (2016)}
\description{
Phylogenetic tree used to generate data matrices in \code{\link{congreveLamsdellMatrices}}.
See Congreve & Lamsdell (2016).
}
\usage{
data(referenceTree)
}
\format{
  A tree of class \code{phylo}
}
\details{

}
\source{
Original source, Sigwart and Sutton (2007); updated by Sutton et al. (2012); modified per Smith & Brazeau (201X)
}
\references{
\itemize{
\item Congreve, C. R., & Lamsdell, J. C. (2016). \cite{Implied weighting and its utility in 
  palaeontological datasets: a study using modelled phylogenetic matrices.} Palaeontology, 
  online ahead of print. doi:10.1111/pala.12236

\item Congreve, C. R., & Lamsdell, J. C. (2016). Data from: Implied weighting and its utility in 
palaeontological datasets: a study using modelled phylogenetic matrices. Dryad Digital
Repository, \href{http://dx.doi.org/10.5061/dryad.7dq0j}{doi:10.5061/dryad.7dq0j}Sigwart, J. D., & Sutton, M. D. (2007). \cite{Deep molluscan phylogeny: synthesis of palaeontological and neontological data.} Proceedings of the Royal Society B: Biological Sciences, 274, 2413-2419. doi:\href{http://dx.doi.org/10.1098/rspb.2007.0701}{10.1098/rspb.2007.0701}
}
}
\examples{
data(referenceTree)
plot(referenceTree)
FitchInfoFast(referenceTree, congreveLamsdellMatrices[[17]])
}
\keyword{datasets}