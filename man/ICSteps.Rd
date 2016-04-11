\name{ICSteps}
\alias{ICSteps}
\title{Information Content Steps}
\description{
  This function estimates the information content of a character \code{char} when e extra steps
  are present, for all possible values of e.
}
\usage{
ICSteps(char, ambiguous.token = 0, expected.minima = 25, max.iter = 10000)
}
\arguments{
  \item{char}{The character in question;}
  \item{ambiguous.token}{Which token, if any, corresponds to the ambigious token (?)
                         (not yet fully implemented);}
  \item{expected.minima}{sample enough trees that the rarest step counts is expected to be seen
                         at least this many times;}
  \item{max.iter}{Maximum iterations to conduct.}
}
\details{
Calculates the number of trees consistent with the character having \emph{e} extra steps, where
e ranges from its minimum possible value (i.e. number of different tokens minus one) to its
maximum.  The number of trees with no extra steps can be calculated exactly; the number of trees
with more additional steps must be approximated.
The function samples \code{n.iter} trees, or enough trees that the trees with the minimum number 
of steps will be recovered at least \code{expected.minima} times, in order to obtain precise 
results.
}

\author{
Martin R. Smith
}
\references{

Faith, D. P. & Trueman, J. W. H. (2001). \cite{Towards an inclusive philosophy for phylogenetic
inference.} Systematic Biology 50:3, 331-350, doi:
\href{http://dx.doi.org/10.1080/10635150118627}{10.1080/10635150118627}

}
\keyword{ tree }

\examples{
# A character that is present in ten taxa and absent in five
character <- c(rep(1, 10), rep(2, 5))
ICSteps (character)
}
