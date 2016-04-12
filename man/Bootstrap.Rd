\name{Bootstrap}
\alias{Bootstrap}
\title{Run a parsimony search on a bootstrapped dataset}
\description{Bootstraps the characters from a dataset and runs a parsimony search.  
Used to explore treespace in the Parsimony Ratchet algorithm (Nixon, 1999)}
\usage{
TreeSearch(tree, data, ParsimonyScorer = ProfileScore, outgroup, concavity = NULL, method = "NNI", maxiter = 100, 
  maxhits = 20, forest.size = 1, cluster = NULL, trace = 1, ...)
}
\arguments{
  \item{tree}{a fully-resolved starting tree in \code{\link{phylo}} format, with the desired outgroup; edge lengths are not supported and will be deleted;}
  \item{data}{a matrix, ideally in the format produced by \code{\link{prepare.data}}, 
  alternatively as a \code{\link{phyDat}} object.  May contain inapplicable data;}
  \item{ParsimonyScorer}{function to use to generate a parsimony score, for example
    \code{\link{parsimony}} or \code{\link{ProfileScore}};}
  \item{outgroup}{a vector listing the taxa in the outgroup;}
  \item{concavity}{concavity constant for implied weighting (not currently implemented!); 
    see \code{\link{parsimony.inapp}};}
  \item{method}{rearrangements to perform; one of \kbd{NNI}, \kbd{SPR}, or \kbd{TBR};}
  \item{maxiter}{the maximum number of iterations to perform before abandoning the search;}
  \item{maxhits}{the maximum times to hit the best pscore before abandoning the search;}
  \item{forest.size}{the maximum number of trees to return - useful in concert with \code{\link{consensus}};}
  \item{cluster}{a cluster prepared using \code{\link{prepare.cluster}}; may speed up search on multicore machines;}
  \item{trace}{higher values provide more verbose user feedback in stdout;}
  \item{\dots}{other arguments to pass to subsequent functions.}
}

\value{
This function returns a tree of class \code{phylo}.  Attributes 'score' and 'hits' will be 
cleared, as the score is a function of the subsampled data and is unrelated to the score 
produced under the original \code{data}.
}
\details{
  Modified from the function in \code{phangorn}, with performance improvements.
}
\references{
\itemize{
 \item Nixon, K. C. (1999). \cite{The Parsimony Ratchet, a new method for rapid parsimony analysis.} 
  Cladistics, 15(4), 407-414. doi:\href{http://dx.doi.org/10.1111/j.1096-0031.1999.tb00277.x}{10.1111/j.1096-0031.1999.tb00277.x}
  
 \item Schliep, K. P. (2011). phangorn: phylogenetic analysis in R. Bioinformatics, 27(4), 592-593.
  doi:\href{http://dx.doi.org/10.1093/bioinformatics/btq706}{10.1093/bioinformatics/btq706}
}
}

\author{
Martin R. Smith
}
\seealso{
\itemize{
\item \code{\link{ProfileScore}}, calculates parsimony score;
\item \code{\link{NNI}}, conducts tree rearrangements;
\item \code{\link{TreeSearch}}, used to search for optimal trees;
\item \code{\link{Pratchet}}, alternative heuristic, useful to escape local optima.
}}
\examples{


\keyword{ tree }
