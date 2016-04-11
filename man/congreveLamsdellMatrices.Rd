\name{SigSut}
\alias{SigSut}
\alias{SigSut.data}
\alias{SigSut.phy}
\alias{SigSut.preparedata}
\alias{SigSut.parstree}
\alias{SigSut.besttree}
\docType{data}
\title{Data from Congreve & Lamsdell 2016}
\description{Contains the phylogenetic matrices from Congreve & Lamsdell (2016)
}
\usage{data(congreveLamsdellMatrices)
}
\format{
  A list of 100 phylogenetic matrices of class phyDat.
}
\details{

}
\source{
Congreve & Lamsdell (2016)
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
data('congreveLamsdellMatrices')
FitchInfoFast(referenceTree, congreveLamsdellMatrices[[17]])
}
\keyword{datasets}