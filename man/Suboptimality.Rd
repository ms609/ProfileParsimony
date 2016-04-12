\name{Suboptimality}
\alias{Suboptimality}

\title{Tree suboptimality extent}
\description{Reports how suboptimal each tree is from the most optimal tree in a set.}
\usage{
Suboptimality(trees)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trees}{trees, of class \code{\link{multiPhylo}}.}
}
\details{
Uses the 'score' attribute of trees generated through the \code{\link{TreeSearch}} function 
(perhaps via \code{\link{Pratchet}}) to list how suboptimal each tree in a set is.
}
\value{A vector of format \code{double(length(trees))}}
\author{
Martin R. Smith
}

\seealso{
\code{\link{TreeSearch}}
\code{\link{Pratchet}}
}
\examples{
data('treeSearchResults')
Suboptimality(treeSearchResults[[1]])
}