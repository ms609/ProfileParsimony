\name{TBR}
\alias{NNI}
\alias{SPR}
\alias{TBR}

\title{Tree rearrangement functions}
\description{These functions performs a single random \acronym{TBR}, \acronym{SPR} or \acronym{NNI} iteration.}
\usage{
NNI(tree)
SPR(tree)
TBR(tree, edge.to.break = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tree}{a fully resolved tree in \code{\link{phyDat}} format;}
  \item{edge.to.break}{the index of an edge to bisect, generated randomly if not specified.}
}
\details{
Performs a single iteration of the nearest-neigbour interchange, subtree pruning and regrafting,
or tree bisection and reconnection algorithms.

NNI and SPR are based on the corresponding phangorn functions, but have been re-coded to 
improve their speed.

Branch lengths are supported by NNI, but not (yet) SPR or TBR.
}
\value{Returns a tree with class \code{phylo}.}
\references{
The algorithms are summarized in

Felsenstein, J. 2004. \cite{Inferring Phylogenies.} Sinauer Associates, Sunderland, Massachusetts.
}
\author{
Martin R. Smith
}

\seealso{
\code{\link{NNI}}
\code{\link{SPR}}
}
\examples{
tree <- rtree(20, br=NULL)
NNI(tree)
SPR(tree)
TBR(tree)
}