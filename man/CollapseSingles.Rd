\name{CollapseSingles}
\alias{CollapseSingles}
\title{Collapse single nodes}
\description{Unknown}
\usage{
Unknown()
}
\arguments{
  \item{arg1}{A phylogenetic tree of class \code{\link{phylo}};}
  \item{arg2}{the node or tip that should form the sister taxon to the new node.  To add a new tip at the root, use "where = 0";}
  \item{arg3}{a character string providing the label the new tip.}
}
\details{
Extends \code{\link{BindTree}}, which cannot handle single-taxon trees.
}
\value{
This function returns a tree of class \code{phylo} with an additional tip at the desired location.
}
\author{Martin R. Smith}

\seealso{
\itemize{
\item \code{\link{BindTree}}
\item \code{\link{nodelabels}}
}
}
\examples{
  plot(tree <- rtree(10, br=NULL)); nodelabels(); nodelabels(15, 15, bg='green'); dev.new()
  plot(AddTip(tree, 15, 'NEW_TIP'))
}
\keyword{ tree }
