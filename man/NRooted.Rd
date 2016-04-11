\name{NRooted}
\alias{NRooted}
\alias{NUnrooted}
\alias{NUnrootedMult}
\alias{LnRooted}
\alias{LnUnrooted}
\alias{LnUnrootedMult}
\title{Number of rooted/unrooted trees}
\description{
  These functions return the number of rooted or unrooted trees consistent with a given pattern
  of splits.
}
\usage{
NRooted(tips)
NUnrooted(tips)
LnRooted(tips)
LnUnrooted(tips)
LnUnrootedMult(splits)
NUnrootedMult(splits)
}
\arguments{
  \item{tips}{The number of tips;}
  \item{splits}{vector listing the number of taxa in each tree bipartition.}
}
\details{

Functions starting N return the number of rooted or unrooted trees, functions starting Ln
provide the log of this number.  Calculations follow Carter et al. 1990, Theorem 2.

}

\author{
Martin R. Smith
}
\references{

Carter, M., Hendy, M., & Penny, D. (1990). \cite{On the distribution of lengths of evolutionary 
trees.} SIAM Journal on Discrete Mathematics, 3(1), 38–47. doi:
\href{http://doi.org/10.1137/0403005}{10.1137/0403005}

}
\keyword{ tree }

\examples{
NRooted(10)
NUnrooted(10)
LnRooted(10)
LnUnrooted(10)
# Number of trees consistent with a character whose states are 00000 11111 222
NUnrootedMult(c(5,5,3))
}
